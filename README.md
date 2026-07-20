# Análisis de Diversidad Funcional

Script en R para calcular índices de diversidad funcional (FRic, FEve, FDiv, FDis) a partir de tablas de muestreos de flora y fauna.

## Requisitos

- R >= 4.1
- RStudio (recomendado) o cualquier entorno R
- Los paquetes se instalan automáticamente la primera vez

## Estructura del repositorio

```
raiz/
├── Proyecto_Tesis.Rproj
├── Tablas_muestreos.ods        ← archivo de datos (4 hojas)
├── Analisis_R/
│   └── analisis_diversidad_funcional_v2.R
└── .gitignore
```

## Uso

1. Clonar el repositorio o descargar los archivos
2. Abrir `Proyecto_Tesis.Rproj` en RStudio
3. Ejecutar `Analisis_R/analisis_diversidad_funcional_v2.R`
4. Los resultados se generan en `Analisis_R/{Flora,Fauna,Unificado}/`

## Formato del archivo de datos (ODS)

El archivo `Tablas_muestreos.ods` debe contener 4 hojas. Los nombres de hojas y columnas son configurables al inicio del script.

### Hojas de abundancia (Flora / Fauna)

| Columna requerida | Descripción |
|---|---|
| Nombre de especie | Se auto-detecta (cualquier nombre con "nombre", "especie", etc.) |
| Zona 1, Zona 2, ... | Columnas numéricas con conteos por zona |
| inventario *(opcional)* | Se excluye del análisis |

Cada fila es un avistamiento/espécimen. Las filas se agrupan por especie y se suman las abundancias por zona.

### Hojas de rasgos

**Flora:**

| Columna | Descripción | Ejemplo |
|---|---|---|
| Nombre_cientifico | Nombre de la especie | *Ceiba pentandra* |
| Altura (m) | Numérico o rango | 15-25 |
| Forma de vida | Categoría | arbol, arbusto, palma, cactus, trepadora |
| Tipo de hoja | Categoría | simple, compuesta, palmada |
| Modo de dispersión | Categoría | viento, animal, gravedad |

**Fauna:**

| Columna | Descripción | Ejemplo |
|---|---|---|
| Nombre_cientifico | Nombre de la especie | *Turdus grayi* |
| Masa corporal (g) | Numérico o rango | 60-80 |
| Dieta principal | Texto | "principalmente insectos" |
| Gremio trófico | Categoría | insectivoro, frugivoro, granivoro |
| Estrato de actividad | Categoría | terrestre, arboreo, volante |

## Configuración para otros datos

Para usar con tablas propias, solo modificar el bloque `CONFIG` al inicio del script:

```r
CONFIG <- list(
  archivo_entrada = "mi_tabla.ods",
  hojas = list(
    flora_abundancia  = "Mi_hoja_flora",
    fauna_abundancia  = "Mi_hoja_fauna",
    rasgos_flora      = "Mis_rasgos_plantas",
    rasgos_fauna      = "Mis_rasgos_animales"
  ),
  columnas = list(
    abundancia = list(
      especie = "Nombre especie",
      excluir = c("observador", "fecha")
    ),
    rasgos_flora = list(
      altura     = "Altura_m",
      forma_vida = "FormaDeVida",
      tipo_hoja  = "Hoja",
      dispersion = "Dispersión"
    ),
    rasgos_fauna = list(
      masa_corporal = "Peso_g",
      dieta         = "Alimentación",
      gremio        = "Gremio",
      estrato       = "Estrato"
    )
  )
)
```

El script valida la estructura del archivo antes de procesar y muestra errores claros si falta alguna columna o hoja.

## Resultados generados

- `indices_{grupo}.csv` — Índices FRic, FEve, FDiv, FDis por zona
- `abundancia_{grupo}.csv` — Matriz de abundancia utilizada
- `pcoa_{grupo}.png` — Gráfico PCoA con zonas coloreadas
- `session_info.txt` — Información de la sesión R para trazabilidad

## Índices calculados

| Índice | Descripción |
|---|---|
| FRic | Riqueza funcional (volumen del convex hull) |
| FEve | Evenness funcional (regularidad del uso del espacio funcional) |
| FDiv | Divergencia funcional (cuánta abundancia se concentra en los rasgos extremos) |
| FDis | Dispersión funcional (distancia media al centroide, ponderada por abundancia) |

---

## Autoría

**Miguel Ángel Morán Piedrahíta**
Biología — Escuela Superior Politécnica del Litoral (ESPOL), Guayaquil, Ecuador
Trabajo de titulación, 2026
mamoran@espol.edu.ec · [linkedin.com/in/miguel-moran-50482b137](https://www.linkedin.com/in/miguel-moran-50482b137)
