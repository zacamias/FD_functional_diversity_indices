###############################################################################
#  ANÁLISIS DE DIVERSIDAD FUNCIONAL v2.2  (REPLICABLE + REUTILIZABLE)
#  Índices: FRic, FEve, FDiv, FDis  (librería FD)
#  Grupos:  Flora, Fauna, Unificado
#  Control de errores estricto — ejecución autónoma de principio a fin
#
#  REPLICABILIDAD:
#    - Rutas relativas a la raíz del repositorio (paquete `here`).
#    - Semilla fija (set.seed) para procesos estocásticos internos.
#    - Impresión de sessionInfo() al final para trazabilidad.
#    - Paquetes instalados solo si faltan.
#
#  REUTILIZACIÓN:
#    - Toda la configuración de entrada está en el bloque CONFIG (sección 1).
#    - Para usar con otros datos, solo modificar las variables del bloque 1.
#    - El script valida la estructura de entrada ANTES de procesar.
#
#  FORMATO ESPERADO DEL ARCHIVO ODS:
#    El archivo debe contener 4 hojas (nombres configurables abajo):
#
#    HOJA DE ABUNDANCIA (Flora / Fauna):
#      - Columna de nombre de especie (cualquier nombre, se auto-detecta).
#      - Columnas de zona: "Zona 1", "Zona 2", ... (nombres con "zona").
#      - Opcional: columna "inventario" (se excluye del análisis).
#      - Cada fila = un espécimen/avistamiento; las filas se agrupan por
#        especie y se suman las abundancias por zona.
#
#    HOJA DE RASGOS FLORA:
#      - Nombre_cientifico (o equivalente, se auto-detecta).
#      - Altura (m): numérico o rango "2-5".
#      - Forma de vida: arbol, arbusto, palma, cactus, trepadora, etc.
#      - Tipo de hoja: simple, compuesta, palmada, etc.
#      - Modo de dispersión: viento, animal, gravedad, etc.
#
#    HOJA DE RASGOS FAUNA:
#      - Nombre_cientifico (o equivalente, se auto-detecta).
#      - Masa corporal (g): numérico o rango.
#      - Dieta principal: "principalmente semillas", "insectos", etc.
#      - Gremio trófico: granivoro, frugivoro, insectivoro, etc.
#      - Estrato de actividad: terrestre, arboreo, volante, etc.
#
#  PARA REPRODUCIR:
#    1. Clonar el repositorio.
#    2. Instalar R >= 4.1.
#    3. Ejecutar este script desde R o RStudio.
#    4. La primera vez se instalan los paquetes faltantes.
#    5. Los resultados se guardan en Analisis_R/{Flora,Fauna,Unificado}/.
###############################################################################

# ── REPRODUCIBILIDAD: SEMILLA ────────────────────────────────────────────────
set.seed(2024)

# ESTRATEGIA ANTE DISTANCIAS NO EUCLÍDEAS:
#   dbFD() requiere que la matriz de distancia sea Euclídea o pasar
#   rasgos originales (él mismo aplica Gower + corrección interna).
#   En este script NO precalculamos la distancia, sino que pasamos
#   la matriz de rasgos directamente a dbFD() con corr="cailliez".
#   Para el análisis Unificado (rasgos con NA cruzados) calculamos
#   una distancia Gower diagonal por bloques, luego extraemos ejes
#   de un PCoA manual (cmdscale) y esos ejes se usan como rasgos
#   sintéticos Euclídeos.
###############################################################################

# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  1. CONFIGURACIÓN — EDITAR AQUÍ PARA OTROS DATOS                       ║
# ╚═══════════════════════════════════════════════════════════════════════════╝
# Solo necesitas modificar esta sección para adaptar el script a tu propia
# tabla de muestreos. Todos los nombres son case-insensitive.

CONFIG <- list(

  # Archivo de entrada (relativo a la raíz del repositorio)
  archivo_entrada = "Tablas_muestreos.ods",

  # Nombres de las hojas en el archivo ODS
  hojas = list(
    flora_abundancia  = "Flora_zonas",
    fauna_abundancia  = "Fauna_zonas",
    rasgos_flora      = "Rasgos_flora",
    rasgos_fauna      = "Rasgos_fauna"
  ),

  # Mapeo de columnas para cada hoja
  # Si tus columnas tienen nombres diferentes, cámbialos aquí.
  # Si el nombre coincide con el esperado, no es necesario cambiarlo.
  columnas = list(
    # Hojas de abundancia: la columna de especie y las columnas a excluir
    abundancia = list(
      especie    = NULL,   # NULL = auto-detectar (busca "nombre", "especie", etc.)
      excluir    = c("inventario")  # columnas a excluir del análisis de zonas
    ),
    # Hojas de rasgos: mapeo nombre_interno -> nombre_en_tu_tabla
    rasgos_flora = list(
      altura         = "Altura (m)",
      forma_vida     = "Forma de vida",
      tipo_hoja      = "Tipo de hoja",
      dispersion     = "Modo de dispersión"
    ),
    rasgos_fauna = list(
      masa_corporal  = "Masa corporal (g)",
      dieta          = "Dieta principal",
      gremio         = "Gremio trófico",
      estrato        = "Estrato de actividad"
    )
  ),

  # Grupos a analizar (comentar/incluir según necesidad)
  analizar = c("Flora", "Fauna", "Unificado")
)

# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  FIN DE LA CONFIGURACIÓN — No editar nada más abajo (a menos que sepas) ║
# ╚═══════════════════════════════════════════════════════════════════════════╝

# ── 2. PAQUETES ──────────────────────────────────────────────────────────────
options(repos = c(CRAN = "https://cloud.r-project.org"))
required_pkgs <- c("here", "FD", "readODS", "tidyverse", "ggrepel",
                    "cluster", "vegan")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
suppressPackageStartupMessages({
  library(here)
  library(FD)
  library(vegan)
  library(readODS)
  library(tidyverse)
  library(ggrepel)
  library(cluster)
})

# ── 3. RUTAS ────────────────────────────────────────────────────────────────
base_dir   <- here()
input_file <- file.path(base_dir, CONFIG$archivo_entrada)
out_root   <- file.path(base_dir, "Analisis_R")

cat("Raíz del proyecto:", base_dir, "\n")
cat("Archivo entrada:  ", input_file, "\n\n")

if (!file.exists(input_file)) {
  stop("No se encontró el archivo de entrada: ", input_file,
       "\nEstructura esperada del repositorio:\n",
       "  raiz/\n",
       "  ├── Analisis_R/\n",
       "  │   └── analisis_diversidad_funcional_v2.R\n",
       "  ├── ", CONFIG$archivo_entrada, "\n",
       "  └── ...")
}

for (sub in c("Flora", "Fauna", "Unificado")) {
  dir.create(file.path(out_root, sub), showWarnings = FALSE, recursive = TRUE)
}

# ── 4. FUNCIONES AUXILIARES ──────────────────────────────────────────────────

extraer_numero <- function(x) {
  v <- as.character(x)
  v <- str_trim(v)
  v <- str_replace_all(v, "~", "")
  v <- str_replace_all(v, ",", "")
  ifelse(str_detect(v, "-"),
    sapply(str_split(v, "-"), function(p) {
      nums <- suppressWarnings(as.numeric(str_trim(p)))
      nums <- nums[!is.na(nums)]
      if (length(nums) >= 2) mean(nums[1:2]) else if (length(nums) == 1) nums[1] else NA_real_
    }),
    suppressWarnings(as.numeric(v))
  )
}

categoria_principal <- function(x) {
  v <- as.character(x)
  v <- str_trim(v)
  v <- str_replace_all(v, "\\s*\\(.*\\)", "")
  v <- str_replace_all(v, "/.*$", "")
  str_trim(v)
}

extraer_fenologia <- function(x) {
  v <- as.character(x)
  case_when(
    str_detect(v, "(?i)caducifolia") ~ "Caducifolia",
    TRUE ~ "Perenne"
  )
}

extraer_tipo_hoja <- function(x) {
  v <- as.character(x)
  tipo <- ifelse(str_detect(v, "/"),
    str_trim(str_extract(v, "[^/]+$")),
    str_trim(v)
  )
  tolower(tipo)
}

limpiar_nombre <- function(df) {
  if (!"Nombre_cientifico" %in% names(df)) {
    idx <- grep("(?i)nombre.*cient|cient|especie|identificacion", names(df))
    if (length(idx) > 0) {
      names(df)[idx[1]] <- "Nombre_cientifico"
    } else {
      stop("No se encontró columna de nombre científico en: ",
           paste(names(df), collapse = ", "))
    }
  }
  df %>%
    mutate(
      Nombre_cientifico = str_squish(as.character(Nombre_cientifico)),
      Nombre_cientifico = ifelse(
        str_detect(Nombre_cientifico, "[a-z]"),
        paste0(toupper(substr(Nombre_cientifico, 1, 1)),
               tolower(substr(Nombre_cientifico, 2, nchar(Nombre_cientifico)))),
        Nombre_cientifico
      )
    ) %>%
    filter(
      !is.na(Nombre_cientifico),
      str_length(Nombre_cientifico) > 0,
      !str_detect(Nombre_cientifico, "(?i)^(indet|sp|sp\\.|spp|spp\\.|sin.*identif|na)$")
    )
}

simplificar_zona <- function(x) {
  str_extract(x, "Zona\\s*\\d+")
}

# ── 4b. VALIDACIÓN DE ESTRUCTURA ─────────────────────────────────────────────
# Función que verifica que el ODS tiene la estructura mínima requerida
# y reporta errores claros al usuario.

validar_estructura <- function(file, cfg) {
  cat("\n--- Validando estructura del archivo de entrada ---\n")
  hojas_disp <- list_ods_sheets(file)
  cat("  Hojas disponibles:", paste(hojas_disp, collapse = ", "), "\n")

  errores <- character()

  # Verificar que las hojas existen
  for (tipo in names(cfg$hojas)) {
    nombre_hoja <- cfg$hojas[[tipo]]
    if (!nombre_hoja %in% hojas_disp) {
      errores <- c(errores,
        paste0("  ✗ Hoja '", nombre_hoja, "' no encontrada (", tipo, "). ",
               "Disponibles: ", paste(hojas_disp, collapse = ", ")))
    }
  }

  if (length(errores) > 0) {
    stop("Errores en la estructura del archivo:\n",
         paste(errores, collapse = "\n"),
         "\n\nRevise los nombres de hojas en CONFIG$hojas.")
  }

  # Verificar columnas en hojas de abundancia
  for (tipo_grupo in c("flora_abundancia", "fauna_abundancia")) {
    nombre_hoja <- cfg$hojas[[tipo_grupo]]
    df <- read_ods(file, sheet = nombre_hoja, limit = 3)

    # Verificar columna de especie
    col_esp <- cfg$columnas$abundancia$especie
    if (!is.null(col_esp)) {
      if (!col_esp %in% names(df)) {
        errores <- c(errores,
          paste0("  ✗ En hoja '", nombre_hoja, "': columna de especie '",
                 col_esp, "' no encontrada. Columnas: ",
                 paste(names(df), collapse = ", ")))
      }
    } else {
      # Auto-detección
      idx <- grep("(?i)nombre.*cient|cient|especie|identificacion", names(df))
      if (length(idx) == 0) {
        errores <- c(errores,
          paste0("  ✗ En hoja '", nombre_hoja,
                 "': no se encontró columna de especie. ",
                 "Use CONFIG$columnas$abundancia$especie para especificarla. ",
                 "Columnas: ", paste(names(df), collapse = ", ")))
      }
    }

    # Verificar que hay columnas de zona
    cols_excluir <- cfg$columnas$abundancia$excluir
    cols_zona <- setdiff(names(df), c("Nombre_cientifico", cols_excluir))
    if (length(cols_zona) == 0) {
      errores <- c(errores,
        paste0("  ✗ En hoja '", nombre_hoja,
               "': no se encontraron columnas de zona. ",
               "Columnas: ", paste(names(df), collapse = ", ")))
    }
  }

  # Verificar columnas en hojas de rasgos
  for (tipo_rasgo in c("rasgos_flora", "rasgos_fauna")) {
    nombre_hoja <- cfg$hojas[[tipo_rasgo]]
    df <- read_ods(file, sheet = nombre_hoja, limit = 3)
    cols_esperadas <- cfg$columnas[[tipo_rasgo]]

    # Verificar columna de especie
    if (!"Nombre_cientifico" %in% names(df)) {
      idx <- grep("(?i)nombre.*cient|cient|especie|identificacion", names(df))
      if (length(idx) == 0) {
        errores <- c(errores,
          paste0("  ✗ En hoja '", nombre_hoja,
                 "': no se encontró columna de nombre científico. ",
                 "Columnas: ", paste(names(df), collapse = ", ")))
      }
    }

    # Verificar rasgos configurados
    for (rasgo in names(cols_esperadas)) {
      col_name <- cols_esperadas[[rasgo]]
      if (!col_name %in% names(df)) {
        # Buscar coincidencia parcial
        coincidencia <- grep("(?i)" %>% paste0(str_split(col_name, " ")[[1]], collapse = "|"),
                             names(df), value = TRUE)
        if (length(coincidencia) == 0) {
          errores <- c(errores,
            paste0("  ✗ En hoja '", nombre_hoja,
                   "': columna de rasgo '", col_name, "' (", rasgo,
                   ") no encontrada. Columnas: ", paste(names(df), collapse = ", ")))
        } else {
          message("  ℹ En hoja '", nombre_hoja, "': columna '", col_name,
                  "' no encontrada. ¿Quiso decir '", coincidencia[1], "'?")
        }
      }
    }
  }

  if (length(errores) > 0) {
    stop("Errores en la estructura del archivo:\n",
         paste(errores, collapse = "\n"))
  }

  cat("  ✓ Estructura válida\n")
}

# ── 5. LECTURA Y VALIDACIÓN DE DATOS ────────────────────────────────────────
cat("\n========================================\n")
cat("  CARGANDO DATOS DESDE ODS\n")
cat("========================================\n\n")

validar_estructura(input_file, CONFIG)

leer_abundancia <- function(file, sheet, cfg) {
  df <- read_ods(file, sheet = sheet)
  # Renombrar columna de especie si está configurada
  col_esp <- cfg$columnas$abundancia$especie
  if (!is.null(col_esp) && col_esp %in% names(df) && col_esp != "Nombre_cientifico") {
    names(df)[names(df) == col_esp] <- "Nombre_cientifico"
  }
  # Excluir columnas configuradas
  cols_excluir <- cfg$columnas$abundancia$excluir
  df <- df %>% select(-any_of(cols_excluir))
  limpiar_nombre(df)
}

flora_abund_raw  <- leer_abundancia(input_file, CONFIG$hojas$flora_abundancia, CONFIG)
fauna_abund_raw  <- leer_abundancia(input_file, CONFIG$hojas$fauna_abundancia, CONFIG)

cat("  Flora abundancia:", nrow(flora_abund_raw), "filas\n")
cat("  Fauna abundancia:", nrow(fauna_abund_raw), "filas\n")

# ── 6. PROCESAR ABUNDANCIAS ──────────────────────────────────────────────────

procesar_abundancia <- function(df, nombre_grupo) {
  cols_zona <- setdiff(names(df), "Nombre_cientifico")
  if (length(cols_zona) == 0) stop("Sin columnas de zona en ", nombre_grupo)
  cat("    Zonas:", paste(cols_zona, collapse = " | "), "\n")

  df %>%
    mutate(across(all_of(cols_zona), ~ suppressWarnings(as.numeric(as.character(.x))))) %>%
    mutate(across(all_of(cols_zona), ~ replace_na(.x, 0))) %>%
    group_by(Nombre_cientifico) %>%
    summarise(across(all_of(cols_zona), ~ sum(.x, na.rm = TRUE), .names = "{.col}"), .groups = "drop")
}

flora_abund <- procesar_abundancia(flora_abund_raw, "Flora")
fauna_abund <- procesar_abundancia(fauna_abund_raw, "Fauna")

cat("  Flora:", nrow(flora_abund), "especies únicas\n")
cat("  Fauna:", nrow(fauna_abund), "especies únicas\n")

# ── 7. PROCESAR RASGOS ───────────────────────────────────────────────────────
# Lee los nombres de columna desde CONFIG y los usa para extraer rasgos.

procesar_rasgos_flora <- function(df, col_cfg) {
  forma_map <- c("arbol" = 1, "arbusto" = 3, "palma" = 2, "cactus" = 4,
                 "trepadora" = 5, "suculenta" = 6, "hierba" = 7, "cicadacea" = 8)

  df %>%
    mutate(
      altura = extraer_numero(.data[[col_cfg$altura]]),
      forma_vida = str_to_lower(str_replace_all(.data[[col_cfg$forma_vida]], "\\s*\\(.*\\)", "")),
      forma_vida = case_when(
        str_detect(forma_vida, "arbusto") ~ "arbusto",
        str_detect(forma_vida, "arbol") ~ "arbol",
        str_detect(forma_vida, "palma") ~ "palma",
        str_detect(forma_vida, "cactus") ~ "cactus",
        str_detect(forma_vida, "trepadora") ~ "trepadora",
        str_detect(forma_vida, "suculenta") ~ "suculenta",
        str_detect(forma_vida, "hierba") ~ "hierba",
        str_detect(forma_vida, "cicad") ~ "cicadacea",
        TRUE ~ "arbol"
      ),
      forma_vida = factor(forma_vida, levels = names(forma_map)),
      fenologia = factor(extraer_fenologia(.data[[col_cfg$tipo_hoja]]),
                         levels = c("Perenne", "Caducifolia")),
      tipo_hoja = extraer_tipo_hoja(.data[[col_cfg$tipo_hoja]]),
      tipo_hoja = case_when(
        str_detect(tipo_hoja, "compuesta|pinnada") ~ "compuesta",
        str_detect(tipo_hoja, "simple") ~ "simple",
        str_detect(tipo_hoja, "palmada") ~ "palmada",
        str_detect(tipo_hoja, "suculenta") ~ "suculenta",
        str_detect(tipo_hoja, "lineal") ~ "lineal",
        str_detect(tipo_hoja, "ancha") ~ "ancha",
        TRUE ~ "simple"
      ),
      tipo_hoja = factor(tipo_hoja),
      dispersion = factor(categoria_principal(.data[[col_cfg$dispersion]]))
    ) %>%
    select(Nombre_cientifico, altura, forma_vida, fenologia, tipo_hoja, dispersion)
}

procesar_rasgos_fauna <- function(df, col_cfg) {
  df %>%
    mutate(
      masa_g = extraer_numero(.data[[col_cfg$masa_corporal]]),
      dieta_raw = .data[[col_cfg$dieta]],
      dieta = str_to_lower(str_extract(dieta_raw, "(?<=(principalmente )).*")),
      dieta = case_when(
        str_detect(dieta, "semilla|gran") ~ "granivoro",
        str_detect(dieta, "fruto|frug") ~ "frugivoro",
        str_detect(dieta, "insecto|artropodo") ~ "insectivoro",
        str_detect(dieta, "hoja|herb") ~ "herbivoro",
        str_detect(dieta, "nectar|polen") ~ "nectarivoro",
        str_detect(dieta_raw, "(?i)semilla") ~ "granivoro",
        str_detect(dieta_raw, "(?i)fruto") ~ "frugivoro",
        str_detect(dieta_raw, "(?i)insecto") ~ "insectivoro",
        str_detect(dieta_raw, "(?i)hoja|herb") ~ "herbivoro",
        str_detect(dieta_raw, "(?i)nectar|polen") ~ "nectarivoro",
        TRUE ~ "omnivoro"
      ),
      dieta = factor(dieta),
      gremio = str_to_lower(str_replace_all(.data[[col_cfg$gremio]], "\\s*\\(.*\\)", "")),
      gremio = factor(gremio),
      estrato = str_to_lower(str_extract(.data[[col_cfg$estrato]], "^[^(]+")),
      estrato = str_trim(estrato),
      estrato = case_when(
        str_detect(estrato, "terrestre.*arboreo|arboreo.*terrestre") ~ "terrestre_arboreo",
        str_detect(estrato, "arboreo") ~ "arboreo",
        str_detect(estrato, "terrestre") ~ "terrestre",
        str_detect(estrato, "volante") ~ "aereo",
        str_detect(estrato, "dosel") ~ "arboreo",
        TRUE ~ "terrestre"
      ),
      estrato = factor(estrato)
    ) %>%
    select(Nombre_cientifico, masa_g, dieta, gremio, estrato)
}

cat("\n--- Procesando rasgos ---\n")

leer_rasgos <- function(file, sheet, col_cfg) {
  df <- read_ods(file, sheet = sheet)
  # Renombrar columna de especie si está configurada
  if (!"Nombre_cientifico" %in% names(df)) {
    idx <- grep("(?i)nombre.*cient|cient|especie|identificacion", names(df))
    if (length(idx) > 0) names(df)[idx[1]] <- "Nombre_cientifico"
  }
  limpiar_nombre(df)
}

rasgos_flora_raw <- leer_rasgos(input_file, CONFIG$hojas$rasgos_flora, CONFIG)
rasgos_fauna_raw <- leer_rasgos(input_file, CONFIG$hojas$rasgos_fauna, CONFIG)

cat("  Rasgos Flora:", nrow(rasgos_flora_raw), "especies\n")
cat("  Rasgos Fauna:", nrow(rasgos_fauna_raw), "especies\n")

rasgos_flora <- procesar_rasgos_flora(rasgos_flora_raw, CONFIG$columnas$rasgos_flora)
rasgos_fauna <- procesar_rasgos_fauna(rasgos_fauna_raw, CONFIG$columnas$rasgos_fauna)

cat("  Flora:", nrow(rasgos_flora), "especies,", ncol(rasgos_flora) - 1, "rasgos\n")
cat("  Fauna:", nrow(rasgos_fauna), "especies,", ncol(rasgos_fauna) - 1, "rasgos\n")

# ── 8. UNIFICAR ABUNDANCIAS ──────────────────────────────────────────────────
cat("\n--- Unificando zonas entre flora y fauna ---\n")

flora_abund_ren <- flora_abund %>%
  rename_with(~ ifelse(str_detect(.x, "(?i)zona"), simplificar_zona(.x), .x))

fauna_abund_ren <- fauna_abund %>%
  rename_with(~ ifelse(str_detect(.x, "(?i)zona"), simplificar_zona(.x), .x))

zonas_flora <- grep("(?i)zona", names(flora_abund_ren), value = TRUE)
zonas_fauna <- grep("(?i)zona", names(fauna_abund_ren), value = TRUE)

all_zonas <- unique(c(zonas_flora, zonas_fauna))
zonas_comunes <- intersect(zonas_flora, zonas_fauna)

if (length(zonas_comunes) == 0 && length(all_zonas) > 0) {
  num_fl <- as.numeric(str_extract(zonas_flora, "\\d+"))
  num_fa <- as.numeric(str_extract(zonas_fauna, "\\d+"))
  num_comunes <- intersect(num_fl, num_fa)
  zonas_comunes <- paste0("Zona ", num_comunes)
}

cat("  Zonas comunes:", paste(zonas_comunes, collapse = ", "), "\n")

if (length(zonas_comunes) == 0) {
  stop("No hay zonas comunes entre flora y fauna para el análisis unificado.")
}

unificado_abund <- bind_rows(
  flora_abund_ren %>% select(Nombre_cientifico, any_of(zonas_comunes)),
  fauna_abund_ren %>% select(Nombre_cientifico, any_of(zonas_comunes))
) %>%
  group_by(Nombre_cientifico) %>%
  summarise(
    across(all_of(zonas_comunes), ~ sum(.x, na.rm = TRUE), .names = "{.col}"),
    .groups = "drop"
  )

cat("  Unificado:", nrow(unificado_abund), "especies\n")

# ── 9. FUNCIÓN PRINCIPAL: dbFD ───────────────────────────────────────────────

analizar_con_dbFD <- function(abund_df, trait_df, label, out_dir) {
  cat("\n", paste(rep("=", 55), collapse = ""), "\n", sep = "")
  cat("  ANÁLISIS:", label, "\n")
  cat(paste(rep("=", 55), collapse = ""), "\n")

  cols_zona <- setdiff(names(abund_df), "Nombre_cientifico")
  sp_comun <- intersect(abund_df$Nombre_cientifico, trait_df$Nombre_cientifico)

  no_trait <- setdiff(abund_df$Nombre_cientifico, trait_df$Nombre_cientifico)
  if (length(no_trait) > 0) {
    cat("  [!]", length(no_trait), "especies sin rasgos:",
        paste(head(no_trait, 5), collapse = ", "),
        if (length(no_trait) > 5) paste("(+", length(no_trait) - 5, "mas)") else "", "\n")
  }

  if (length(sp_comun) < 3) {
    cat("  [!] < 3 especies con rasgos. SE OMITE.\n")
    return(NULL)
  }

  abund <- abund_df %>%
    filter(Nombre_cientifico %in% sp_comun) %>%
    arrange(Nombre_cientifico)

  abund_matrix <- abund %>%
    column_to_rownames("Nombre_cientifico") %>%
    as.matrix() %>%
    t()

  trait <- trait_df %>%
    filter(Nombre_cientifico %in% sp_comun) %>%
    arrange(Nombre_cientifico)
  trait_matrix <- trait %>% column_to_rownames("Nombre_cientifico")

  for (i in seq_len(ncol(trait_matrix))) {
    if (is.character(trait_matrix[, i])) trait_matrix[, i] <- as.factor(trait_matrix[, i])
  }

  n_sp <- ncol(abund_matrix)
  n_traits <- ncol(trait_matrix)
  n_sites <- nrow(abund_matrix)

  cat("  Especies:", n_sp, "| Rasgos:", n_traits, "| Zonas:", n_sites, "\n")

  # Filtros
  zero_sp <- colSums(abund_matrix, na.rm = TRUE) == 0
  if (any(zero_sp)) {
    cat("  [!]", sum(zero_sp), "especies sin avistamientos eliminadas\n")
    abund_matrix <- abund_matrix[, !zero_sp, drop = FALSE]
    trait_matrix <- trait_matrix[!zero_sp, , drop = FALSE]
    n_sp <- ncol(abund_matrix)
  }

  zero_zone <- rowSums(abund_matrix, na.rm = TRUE) == 0
  if (any(zero_zone)) {
    cat("  [!]", sum(zero_zone), "zonas sin especies:",
        paste(rownames(abund_matrix)[zero_zone], collapse = ", "), "\n")
    abund_matrix <- abund_matrix[!zero_zone, , drop = FALSE]
    n_sites <- nrow(abund_matrix)
  }

  constant_trait <- apply(trait_matrix, 2, function(x) {
    if (is.numeric(x)) sd(x, na.rm = TRUE) == 0 || all(is.na(x))
    else length(unique(x[!is.na(x)])) <= 1
  })
  if (any(constant_trait)) {
    cat("  [!]", sum(constant_trait), "rasgos constantes eliminados\n")
    trait_matrix <- trait_matrix[, !constant_trait, drop = FALSE]
    n_traits <- ncol(trait_matrix)
  }

  if (n_sp < 3) { cat("  [!] < 3 especies. SE OMITE.\n"); return(NULL) }
  if (n_sites < 1) { cat("  [!] 0 zonas. SE OMITE.\n"); return(NULL) }
  if (n_traits < 1) { cat("  [!] 0 rasgos. SE OMITE.\n"); return(NULL) }

  calc_fric <- n_sp > n_traits + 1
  m_val <- if (calc_fric) n_traits else max(n_sp - 2, 2)
  if (!calc_fric) {
    cat("  [!] n_sp (", n_sp, ") <= n_traits (", n_traits,
        "): calc.FRic = FALSE, m =", m_val, "\n")
  }

  # Fallback con 4 intentos
  ejecutar <- function() {
    tryCatch({
      cat("  -> Intento 1: dbFD con rasgos + corr=cailliez\n")
      dbFD(x = trait_matrix, a = abund_matrix,
           calc.FRic = calc_fric, m = m_val,
           stand.FRic = TRUE, scale.RaoQ = TRUE,
           calc.CWM = FALSE, calc.FGR = FALSE,
           clust.type = "ward", messages = FALSE,
           corr = "cailliez")
    }, error = function(e1) {
      cat("  [!] Intento 1 falló:", e1$message, "\n")
      tryCatch({
        cat("  -> Intento 2: dbFD con rasgos + corr=lingoes\n")
        dbFD(x = trait_matrix, a = abund_matrix,
             calc.FRic = calc_fric, m = m_val,
             stand.FRic = TRUE, scale.RaoQ = TRUE,
             calc.CWM = FALSE, calc.FGR = FALSE,
             clust.type = "ward", messages = FALSE,
             corr = "lingoes")
      }, error = function(e2) {
        cat("  [!] Intento 2 falló:", e2$message, "\n")
        tryCatch({
          cat("  -> Intento 3: PCoA manual + rasgos sintéticos\n")
          dist_g <- gowdis(trait_matrix)
          k_axes <- max(min(n_sp - 2, n_traits * 2, 10), 2)
          pcoa <- cmdscale(dist_g, k = k_axes, eig = TRUE)
          n_cols <- ncol(pcoa$points)
          pos <- which(pcoa$eig[1:n_cols] > 1e-10)
          if (length(pos) < 2) pos <- 1:min(2, n_cols)
          syn_traits <- pcoa$points[, pos, drop = FALSE]
          colnames(syn_traits) <- paste0("PCoA_", seq_len(ncol(syn_traits)))
          n_syn <- ncol(syn_traits)
          cf <- n_sp > n_syn + 1
          m <- if (cf) n_syn else max(n_sp - 2, 2)
          if (!cf) cat("    calc.FRic = FALSE con rasgos sintéticos\n")
          dbFD(x = syn_traits, a = abund_matrix,
               calc.FRic = cf, m = m,
               stand.FRic = TRUE, scale.RaoQ = TRUE,
               calc.CWM = FALSE, calc.FGR = FALSE,
               clust.type = "ward", messages = FALSE)
        }, error = function(e3) {
          cat("  [!] Intento 3 falló:", e3$message, "\n")
          tryCatch({
            cat("  -> Intento 4: PCoA + calc.FRic = FALSE\n")
            dist_g <- gowdis(trait_matrix)
            k_axes <- max(min(n_sp - 2, 5), 2)
            pcoa <- cmdscale(dist_g, k = k_axes, eig = TRUE)
            n_cols <- ncol(pcoa$points)
            pos <- which(pcoa$eig[1:n_cols] > 1e-10)
            if (length(pos) < 2) pos <- 1:min(2, n_cols)
            syn_traits <- pcoa$points[, pos, drop = FALSE]
            colnames(syn_traits) <- paste0("PCoA_", seq_len(ncol(syn_traits)))
            dbFD(x = syn_traits, a = abund_matrix,
                 calc.FRic = FALSE, m = 2,
                 stand.FRic = TRUE, scale.RaoQ = TRUE,
                 calc.CWM = FALSE, calc.FGR = FALSE,
                 clust.type = "ward", messages = FALSE)
          }, error = function(e4) {
            cat("  [!] Todos los intentos fallaron. SE OMITE.\n")
            NULL
          })
        })
      })
    })
  }

  resultado <- ejecutar()
  if (is.null(resultado)) return(NULL)

  # Extraer índices
  indices <- data.frame(
    Zona = rownames(abund_matrix),
    FRic = if (!is.null(resultado$FRic)) resultado$FRic else NA,
    FEve = if (!is.null(resultado$FEve)) resultado$FEve else NA,
    FDiv = if (!is.null(resultado$FDiv)) resultado$FDiv else NA,
    FDis = if (!is.null(resultado$FDis)) resultado$FDis else NA,
    stringsAsFactors = FALSE
  )

  cat("\n  Índices:\n")
  print(indices, row.names = FALSE)

  # Guardar
  csv_path <- file.path(out_dir, paste0("indices_", tolower(label), ".csv"))
  write.csv(indices, csv_path, row.names = FALSE)
  cat("  CSV:", csv_path, "\n")

  abund_csv <- file.path(out_dir, paste0("abundancia_", tolower(label), ".csv"))
  write.csv(abund_matrix, abund_csv, row.names = TRUE)
  cat("  Abundancia:", abund_csv, "\n")

  # Gráfico PCoA
  tryCatch({
    pcoa_coord <- NULL
    var1 <- NA; var2 <- NA

    if (!is.null(resultado$x) && ncol(resultado$x) >= 2) {
      pcoa_coord <- resultado$x[, 1:2]
      var1 <- if (!is.null(resultado$prop_expli) && length(resultado$prop_expli) >= 1)
        round(100 * resultado$prop_expli[1], 1) else NA
      var2 <- if (!is.null(resultado$prop_expli) && length(resultado$prop_expli) >= 2)
        round(100 * resultado$prop_expli[2], 1) else NA
    } else {
      if (nrow(abund_matrix) >= 3) {
        bc_dist <- vegdist(abund_matrix, method = "bray", na.rm = TRUE)
        pcoa_m <- cmdscale(bc_dist, k = 2, eig = TRUE)
        n_cols <- ncol(pcoa_m$points)
        if (n_cols >= 2) {
          pcoa_coord <- pcoa_m$points[, 1:2]
          ve <- abs(pcoa_m$eig) / sum(abs(pcoa_m$eig))
          var1 <- round(100 * ve[1], 1)
          var2 <- round(100 * ve[2], 1)
        }
      }
    }

    if (!is.null(pcoa_coord) && ncol(pcoa_coord) >= 2) {
      pcoa_df <- as.data.frame(pcoa_coord)
      colnames(pcoa_df) <- c("PCoA1", "PCoA2")
      pcoa_df$Zona <- rownames(abund_matrix)

      if (!is.null(resultado$FRic) && length(resultado$FRic) == nrow(pcoa_df)) {
        pcoa_df$FRic <- resultado$FRic
        p <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, label = Zona, color = FRic))
      } else {
        p <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, label = Zona))
      }

      p <- p +
        geom_point(size = 4, alpha = 0.85) +
        geom_text_repel(size = 4.5, force = 3, max.overlaps = 20,
                        box.padding = 0.5, point.padding = 0.3) +
        labs(
          title = paste("PCoA —", label),
          subtitle = "Diversidad funcional (FRic, FEve, FDiv, FDis)",
          x = if (!is.na(var1)) paste0("PCoA1 (", var1, "%)") else "PCoA1",
          y = if (!is.na(var2)) paste0("PCoA2 (", var2, "%)") else "PCoA2"
        ) +
        theme_classic(base_size = 14) +
        theme(
          plot.title    = element_text(face = "bold", size = 16, hjust = 0.5),
          plot.subtitle = element_text(size = 11, hjust = 0.5),
          axis.title    = element_text(size = 13),
          axis.text     = element_text(size = 11),
          legend.position = "right"
        )

      if (!is.null(resultado$FRic)) {
        p <- p + scale_color_gradient(low = "blue", high = "red", na.value = "grey50")
      }

      png_path <- file.path(out_dir, paste0("pcoa_", tolower(label), ".png"))
      ggsave(png_path, plot = p, width = 10, height = 8, dpi = 300)
      cat("  Gráfico:", png_path, "\n")
    } else {
      cat("  PCoA no disponible (< 2 dimensiones). Gráfico omitido.\n")
    }
  }, error = function(e) {
    cat("  [!] Gráfico PCoA falló:", e$message, "- omitido.\n")
  })

  cat(paste(rep("=", 55), collapse = ""), "\n")
  invisible(indices)
}

# ── 10. VERSIÓN UNIFICADO ────────────────────────────────────────────────────

analizar_unificado <- function(abund_df, traits_flora, traits_fauna, label, out_dir) {
  cat("\n", paste(rep("=", 55), collapse = ""), "\n", sep = "")
  cat("  ANÁLISIS:", label, "\n")
  cat(paste(rep("=", 55), collapse = ""), "\n")

  sp_tot <- abund_df$Nombre_cientifico
  sp_fl <- intersect(sp_tot, traits_flora$Nombre_cientifico)
  sp_fa <- intersect(sp_tot, traits_fauna$Nombre_cientifico)

  cat("  Especies flora con rasgos:", length(sp_fl), "\n")
  cat("  Especies fauna con rasgos:", length(sp_fa), "\n")

  if (length(sp_fl) + length(sp_fa) < 3) {
    cat("  [!] < 3 especies con rasgos. SE OMITE.\n")
    return(NULL)
  }

  d_g_fl <- NULL
  if (length(sp_fl) > 1) {
    tf <- traits_flora %>% filter(Nombre_cientifico %in% sp_fl) %>% arrange(Nombre_cientifico)
    tm <- tf %>% column_to_rownames("Nombre_cientifico")
    d_g_fl <- as.matrix(gowdis(tm))
  }

  d_g_fa <- NULL
  if (length(sp_fa) > 1) {
    tf <- traits_fauna %>% filter(Nombre_cientifico %in% sp_fa) %>% arrange(Nombre_cientifico)
    tm <- tf %>% column_to_rownames("Nombre_cientifico")
    d_g_fa <- as.matrix(gowdis(tm))
  }

  spp <- unique(c(sp_fl, sp_fa))
  n <- length(spp)
  mat_comb <- matrix(1, nrow = n, ncol = n, dimnames = list(spp, spp))
  diag(mat_comb) <- 0

  if (!is.null(d_g_fl)) mat_comb[sp_fl, sp_fl] <- d_g_fl[sp_fl, sp_fl]
  if (!is.null(d_g_fa)) mat_comb[sp_fa, sp_fa] <- d_g_fa[sp_fa, sp_fa]

  abund <- abund_df %>%
    filter(Nombre_cientifico %in% spp) %>%
    arrange(Nombre_cientifico)

  abund_matrix <- abund %>%
    column_to_rownames("Nombre_cientifico") %>%
    as.matrix() %>%
    t()

  sp_en_abund <- colnames(abund_matrix)
  sp_en_dist <- rownames(mat_comb)
  sp_final <- intersect(sp_en_abund, sp_en_dist)

  if (length(sp_final) < 3) {
    cat("  [!] < 3 especies tras intersección. SE OMITE.\n")
    return(NULL)
  }

  abund_matrix <- abund_matrix[, sp_final, drop = FALSE]
  mat_comb <- mat_comb[sp_final, sp_final]

  cat("  Especies:", ncol(abund_matrix), "| Zonas:", nrow(abund_matrix), "\n")

  zero_sp <- colSums(abund_matrix, na.rm = TRUE) == 0
  if (any(zero_sp)) {
    cat("  [!]", sum(zero_sp), "especies sin avistamientos eliminadas\n")
    abund_matrix <- abund_matrix[, !zero_sp, drop = FALSE]
    mat_comb <- mat_comb[!zero_sp, !zero_sp]
  }

  zero_zone <- rowSums(abund_matrix, na.rm = TRUE) == 0
  if (any(zero_zone)) {
    cat("  [!]", sum(zero_zone), "zonas sin especies eliminadas\n")
    abund_matrix <- abund_matrix[!zero_zone, , drop = FALSE]
  }

  if (ncol(abund_matrix) < 3) { cat("  [!] < 3 especies. SE OMITE.\n"); return(NULL) }
  if (nrow(abund_matrix) < 1) { cat("  [!] 0 zonas. SE OMITE.\n"); return(NULL) }

  n_sp <- ncol(abund_matrix)
  k_axes <- max(min(n_sp - 2, 8), 2)
  pcoa <- cmdscale(as.dist(mat_comb), k = k_axes, eig = TRUE)
  n_cols <- ncol(pcoa$points)
  pos <- which(pcoa$eig[1:n_cols] > 1e-10)
  if (length(pos) < 2) pos <- 1:min(2, n_cols)
  syn_traits <- pcoa$points[, pos, drop = FALSE]
  colnames(syn_traits) <- paste0("PCoA_", seq_len(ncol(syn_traits)))

  n_syn <- ncol(syn_traits)
  calc_fric <- n_sp > n_syn + 1
  m_val <- if (calc_fric) n_syn else max(n_sp - 2, 2)
  if (!calc_fric) cat("  [!] calc.FRic = FALSE (pocas especies para", n_syn, "ejes)\n")

  resultado <- tryCatch(
    dbFD(x = syn_traits, a = abund_matrix, calc.FRic = calc_fric, m = m_val,
         stand.FRic = TRUE, scale.RaoQ = TRUE, calc.CWM = FALSE,
         calc.FGR = FALSE, clust.type = "ward", messages = FALSE),
    error = function(e) {
      cat("  [!] dbFD falló:", e$message, "\n")
      cat("  -> Reintentando con calc.FRic = FALSE\n")
      tryCatch(
        dbFD(x = syn_traits, a = abund_matrix, calc.FRic = FALSE, m = 2,
             stand.FRic = TRUE, scale.RaoQ = TRUE, calc.CWM = FALSE,
             calc.FGR = FALSE, clust.type = "ward", messages = FALSE),
        error = function(e2) {
          cat("  [!] Unificado falló:", e2$message, "\n")
          NULL
        }
      )
    }
  )

  if (is.null(resultado)) return(NULL)

  indices <- data.frame(
    Zona = rownames(abund_matrix),
    FRic = if (!is.null(resultado$FRic)) resultado$FRic else NA,
    FEve = if (!is.null(resultado$FEve)) resultado$FEve else NA,
    FDiv = if (!is.null(resultado$FDiv)) resultado$FDiv else NA,
    FDis = if (!is.null(resultado$FDis)) resultado$FDis else NA,
    stringsAsFactors = FALSE
  )

  cat("\n  Índices:\n")
  print(indices, row.names = FALSE)

  csv_path <- file.path(out_dir, paste0("indices_", tolower(label), ".csv"))
  write.csv(indices, csv_path, row.names = FALSE)
  cat("  CSV:", csv_path, "\n")

  # Gráfico
  tryCatch({
    pcoa_coord <- NULL
    var1 <- NA; var2 <- NA

    if (!is.null(resultado$x) && ncol(resultado$x) >= 2) {
      pcoa_coord <- resultado$x[, 1:2]
      var1 <- if (!is.null(resultado$prop_expli) && length(resultado$prop_expli) >= 1)
        round(100 * resultado$prop_expli[1], 1) else NA
      var2 <- if (!is.null(resultado$prop_expli) && length(resultado$prop_expli) >= 2)
        round(100 * resultado$prop_expli[2], 1) else NA
    } else {
      if (nrow(abund_matrix) >= 3) {
        bc_dist <- vegdist(abund_matrix, method = "bray", na.rm = TRUE)
        pcoa_m <- cmdscale(bc_dist, k = 2, eig = TRUE)
        n_cols <- ncol(pcoa_m$points)
        if (n_cols >= 2) {
          pcoa_coord <- pcoa_m$points[, 1:2]
          ve <- abs(pcoa_m$eig) / sum(abs(pcoa_m$eig))
          var1 <- round(100 * ve[1], 1)
          var2 <- round(100 * ve[2], 1)
        }
      }
    }

    if (!is.null(pcoa_coord) && ncol(pcoa_coord) >= 2) {
      pcoa_df <- as.data.frame(pcoa_coord)
      colnames(pcoa_df) <- c("PCoA1", "PCoA2")
      pcoa_df$Zona <- rownames(abund_matrix)

      p <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, label = Zona)) +
        geom_point(size = 4, alpha = 0.85, color = "#2166ac") +
        geom_text_repel(size = 4.5, force = 3, max.overlaps = 20,
                        box.padding = 0.5, point.padding = 0.3) +
        labs(
          title = paste("PCoA —", label),
          subtitle = "Diversidad funcional unificada (Flora + Fauna)",
          x = if (!is.na(var1)) paste0("PCoA1 (", var1, "%)") else "PCoA1",
          y = if (!is.na(var2)) paste0("PCoA2 (", var2, "%)") else "PCoA2"
        ) +
        theme_classic(base_size = 14) +
        theme(
          plot.title    = element_text(face = "bold", size = 16, hjust = 0.5),
          plot.subtitle = element_text(size = 11, hjust = 0.5),
          axis.title    = element_text(size = 13),
          axis.text     = element_text(size = 11)
        )

      png_path <- file.path(out_dir, paste0("pcoa_", tolower(label), ".png"))
      ggsave(png_path, plot = p, width = 10, height = 8, dpi = 300)
      cat("  Gráfico:", png_path, "\n")
    }
  }, error = function(e) {
    cat("  [!] Gráfico PCoA falló:", e$message, "- omitido.\n")
  })

  cat(paste(rep("=", 55), collapse = ""), "\n")
  invisible(indices)
}

# ── 11. EJECUTAR ─────────────────────────────────────────────────────────────

resultados <- list()

cat("\n\n========== EJECUTANDO ANÁLISIS ==========\n")

if ("Flora" %in% CONFIG$analizar) {
  resultados[["Flora"]] <- tryCatch(
    analizar_con_dbFD(flora_abund, rasgos_flora, "Flora",
                      file.path(out_root, "Flora")),
    error = function(e) { cat("\n[!] ERROR Flora:", e$message, "\n"); NULL }
  )
}

if ("Fauna" %in% CONFIG$analizar) {
  resultados[["Fauna"]] <- tryCatch(
    analizar_con_dbFD(fauna_abund, rasgos_fauna, "Fauna",
                      file.path(out_root, "Fauna")),
    error = function(e) { cat("\n[!] ERROR Fauna:", e$message, "\n"); NULL }
  )
}

if ("Unificado" %in% CONFIG$analizar) {
  resultados[["Unificado"]] <- tryCatch(
    analizar_unificado(unificado_abund, rasgos_flora, rasgos_fauna,
                       "Unificado", file.path(out_root, "Unificado")),
    error = function(e) { cat("\n[!] ERROR Unificado:", e$message, "\n"); NULL }
  )
}

# ── 12. RESUMEN ─────────────────────────────────────────────────────────────
cat("\n\n========================================\n")
cat("  RESUMEN DE RESULTADOS\n")
cat("========================================\n\n")

for (grp in names(resultados)) {
  cat("  [", grp, "]\n", sep = "")
  if (!is.null(resultados[[grp]])) {
    cat("    Zonas:", nrow(resultados[[grp]]), "\n")
    cat("    CSV:  ", file.path(out_root, grp, paste0("indices_", tolower(grp), ".csv")), "\n")
    png_path <- file.path(out_root, grp, paste0("pcoa_", tolower(grp), ".png"))
    if (file.exists(png_path)) cat("    PNG:  ", png_path, "\n")
  } else {
    cat("    [OMITIDO]\n")
  }
}

cat("\n========================================\n")
cat("  ANÁLISIS COMPLETADO\n")
cat("========================================\n")

# ── 13. EXPORTAR PLANTILLAS Y RASGOS ────────────────────────────────────────

cat("\n--- Exportando plantillas vacías y rasgos completos ---\n")

plantilla_flora <- data.frame(
  Nombre_cientifico = flora_abund$Nombre_cientifico,
  altura_m          = NA_real_,
  forma_vida        = factor(NA, levels = names(c(
    "arbol"=1,"arbusto"=3,"palma"=2,"cactus"=4,
    "trepadora"=5,"suculenta"=6,"hierba"=7,"cicadacea"=8))),
  tipo_hoja         = factor(NA, levels = c("simple","compuesta","palmada",
    "suculenta","lineal","ancha")),
  fenologia         = factor(NA, levels = c("Perenne","Caducifolia")),
  dispersion        = factor(NA, levels = c("anemocora","zoocora","autocora")),
  stringsAsFactors  = FALSE
)

plantilla_fauna <- data.frame(
  Nombre_cientifico = fauna_abund$Nombre_cientifico,
  masa_g            = NA_real_,
  dieta             = factor(NA, levels = c("granivoro","frugivoro",
    "insectivoro","herbivoro","nectarivoro","omnivoro")),
  gremio            = factor(NA),
  estrato           = factor(NA, levels = c("arboreo","terrestre",
    "aereo","terrestre_arboreo")),
  stringsAsFactors  = FALSE
)

plantilla_unificado <- bind_rows(
  plantilla_flora %>% mutate(Grupo = "Flora"),
  plantilla_fauna %>% mutate(Grupo = "Fauna")
) %>% select(Grupo, everything())

cat(sprintf("  Plantilla Flora:     %d especies\n", nrow(plantilla_flora)))
cat(sprintf("  Plantilla Fauna:     %d especies\n", nrow(plantilla_fauna)))
cat(sprintf("  Plantilla Unificado: %d especies\n", nrow(plantilla_unificado)))

write.csv(plantilla_flora,
          file.path(out_root, "Flora", "plantilla_rasgos_flora.csv"),
          row.names = FALSE)
write.csv(plantilla_fauna,
          file.path(out_root, "Fauna", "plantilla_rasgos_fauna.csv"),
          row.names = FALSE)
write.csv(plantilla_unificado,
          file.path(out_root, "Unificado", "plantilla_rasgos_unificado.csv"),
          row.names = FALSE)

write.csv(rasgos_flora,
          file.path(out_root, "Flora", "rasgos_completos_flora.csv"),
          row.names = FALSE)
write.csv(rasgos_fauna,
          file.path(out_root, "Fauna", "rasgos_completos_fauna.csv"),
          row.names = FALSE)

cat("  Plantillas y rasgos exportados a cada subdirectorio.\n")

# ── 14. INFORMACIÓN DE SESIÓN ───────────────────────────────────────────────
cat("\n--- sessionInfo() ---\n")
session_info <- capture.output(sessionInfo())
cat(paste(session_info, collapse = "\n"), "\n")

si_path <- file.path(out_root, "session_info.txt")
writeLines(session_info, si_path)
cat("\nSession info guardada en:", si_path, "\n")
