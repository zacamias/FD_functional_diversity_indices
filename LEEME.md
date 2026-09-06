# Diversidad funcional del campus Las Peñas (ESPOL)

[![Licencia: MIT](https://img.shields.io/badge/Licencia-MIT-yellow.svg)](LICENSE)
[![Datos: CC BY 4.0](https://img.shields.io/badge/Datos-CC%20BY%204.0-lightgrey.svg)](LICENSE)
[![R](https://img.shields.io/badge/R-%E2%89%A5%204.3-276DC3.svg)](https://www.r-project.org/)

Versión en español del [README](README.md).

Pipeline de R y conjunto de datos de campo que caracterizan la **diversidad
funcional** de la flora y fauna del campus Las Peñas (ESPOL, Guayaquil,
Ecuador) como **línea de referencia previa** a la remodelación del campus.

Proyecto Integrador, Biología · Miguel Ángel Morán Piedrahita · ESPOL, 2026.
Tutor: Julián Alfredo Pérez Correa, PhD.

---

## Por qué existe

Más de la mitad de un campus de 3.25 ha —situado entre los últimos fragmentos de
bosque seco dentro de Guayaquil, a pocos metros del río Guayas— está siendo
reconstruido, y nadie había registrado qué vivía allí. El bosque seco estacional
del Ecuador perdió 2 631.91 km² entre 1990 y 2018; apenas el 13 % del remanente
costero tiene alguna figura de protección, y cerca del 80 % de su biota es
endémica del Centro de Endemismo Tumbesino.

Contar especies no responde la pregunta que plantea una obra civil. La riqueza no
dice **qué hacen** esas especies, y es eso lo que la obra interrumpe. Este
repositorio mide los rasgos.

## Qué calcula

Por zona y por grupo taxonómico:

| Índice | Qué captura |
|---|---|
| **FRic** | Riqueza funcional — volumen del espacio de rasgos ocupado |
| **FEve** | Uniformidad funcional — regularidad de las abundancias en ese volumen |
| **FDiv** | Divergencia funcional — cuán lejos del centro están las especies abundantes |
| **FDis** | Dispersión funcional — distancia media al centroide, ponderada por abundancia |
| **Q de Rao** | Entropía cuadrática — distancia media entre pares, ponderada por abundancia |
| **FR** | Redundancia funcional — fracción de la diversidad de especies que *no* es diversidad funcional |
| **SES** | Lo anterior estandarizado contra un modelo nulo, controlando la riqueza |

## Resultados

![Diversidad funcional de la flora por zona](figures/heatmap_flora.png)

![Diversidad funcional de las aves por zona](figures/heatmap_birds.png)

Dos patrones dominan y van juntos: la **riqueza funcional es baja en todas las
zonas** —la compresión que produce de forma sistemática el filtrado urbano,
acentuada por una flora que es 100 % introducida y ornamental— mientras la
**divergencia funcional es alta** (0.59–1.00): las especies abundantes se sitúan
en la periferia de ese volumen reducido.

Un resultado merece lectura cuidadosa. En aves, FRic es `NA` para las Zonas 3, 4
y 5 porque contienen *exactamente las mismas cuatro especies*, y un volumen sin
ponderar no puede distinguir listas idénticas. La Zona 1 es la excepción: la
única con *Columba livia* y *Quiscalus mexicanus*. Que el **estacionamiento**
resulte así la zona de mayor riqueza funcional de aves del campus es el mejor
argumento disponible contra usar un índice aislado como criterio de
conservación.

![La diversidad de especies no es diversidad funcional](figures/taxonomic_vs_functional.png)

La brecha entre las dos barras es la redundancia funcional: especies cuyo papel
ecológico ya está cubierto por otra presente. Es la medición que sostiene el
argumento de resiliencia — un FRic bajo dice que la comunidad ocupa poco espacio
de rasgos, pero **no** dice que distintas especies cumplan la misma función. Esa
es otra afirmación y necesita su propio número.

![Redundancia funcional por zona](figures/redundancy_by_zone.png)

⚠️ **Las dos series de redundancia no coinciden.** Es un punto de verificación
abierto y documentado ([`docs/verification.md`](docs/verification.md), ítem
**V1**), no un resultado cerrado. Ambas usan `FR = 1 − Q/D`; la diferencia es
sistemática y muy probablemente se debe a cómo se escaló la disimilitud dentro de
`Q`. No cite una cifra de redundancia de este repositorio sin leer antes esa
entrada.

## Uso rápido

```bash
git clone https://github.com/zacamias/FD_functional_diversity_indices.git
cd FD_functional_diversity_indices

# Figuras y análisis de redundancia / diversidad taxonómica
# (solo necesita ggplot2, reshape2, cluster, vegan — no FD)
Rscript R/figures_reference.R

# Pipeline completo: todos los índices, modelos nulos, sensibilidad (necesita FD)
Rscript R/functional_diversity_analysis.R
```

Ejecute ambos **desde la raíz del repositorio**, no desde `R/`. Las salidas van a
`results/`, que está en `.gitignore` y se regenera en cada corrida.

Para restaurar las versiones exactas de los paquetes:

```r
install.packages("renv")
renv::restore()
```

> `renv.lock` se redactó a mano a partir de las versiones reportadas en la tesis
> (§2.9). Ejecute `renv::snapshot()` en su propia máquina para reemplazarlo por
> un lockfile con hashes verificados — esa es la versión confiable.

## Método

| Decisión | Elección |
|---|---|
| Abundancia de flora | Cobertura por proyección de copa, fracción 0–1 |
| Abundancia de fauna | Detecciones por zona y campaña |
| Distancia | Gower (rasgos mixtos, reescalado interno por rango) |
| Transformación | `log10` en rasgos continuos sesgados: masa, altura, SLA |
| Ordenación | PCoA, corrección de **Cailliez** para autovalores negativos |
| Motor | `FD::dbFD`, con `mFD` como verificación cruzada opcional |
| Redundancia | Ricotta et al. (2016), `FR = D − Q`, reportada como `1 − Q/D` |
| Modelo nulo | 999 permutaciones de etiquetas de rasgos → SES + *p* exacta bilateral |
| Sensibilidad | Índices recalculados con abundancia cruda, `log(x+1)` y presencia/ausencia |

**Sobre la corrección.** Las distancias de Gower sobre rasgos mixtos
generalmente no son euclidianas, lo que produce autovalores negativos en la PCoA
y deja mal definida la envolvente convexa que sostiene FRic. Lingoes añade una
constante menor y distorsiona menos, así que en el fondo es la mejor opción
técnica. **Aun así, Cailliez es el valor por defecto aquí**, porque es lo que
produjo las tablas publicadas, y un repositorio cuya función es reproducirlas
debe partir de lo que las produjo. Cambiar es una línea, y pertenece a un
análisis de sensibilidad reportado, no a una edición silenciosa.

## Por qué las aves se analizan por separado

*Apis mellifera* representó entre el **90.1 % y el 95.7 %** de los registros de
fauna de la Zona 4. Como FDis pondera las distancias por la abundancia relativa,
integrarla con los vertebrados arrastra el centroide de la comunidad hacia un
único tipo funcional y deprime el índice: un artefacto metodológico, no un
hallazgo ecológico. Al analizar las aves por separado el efecto desaparece.

Dos razones adicionales hacen la separación obligatoria y no meramente
conveniente. Los conteos de la abeja son **extrapolados**, no observados: cuatro
unidades de 2.5 m en jardineras de *Ixora coccinea*, escaladas por un factor de
30.4 a los 304 m que ocupa la planta. Y esas unidades se ubicaron en las cuatro
jardineras *más largas*, no al azar — como las abejas concentran forrajeras en
los parches más ricos mediante reclutamiento por danza, la cifra es un **límite
superior**: un índice de intensidad de visitación por unidad de recurso floral, no
un tamaño poblacional.

## Control de la riqueza

Todos estos índices escalan con la riqueza por construcción, y las zonas van de
2 a 26 especies de flora. Comparar zonas directamente confunde *cuántas especies
hay* con *cómo se distribuyen en el espacio de rasgos*: decir "la Zona 4 es
funcionalmente más rica" corre el riesgo de no significar más que "la Zona 4
tiene más especies".

`R/null_models.R` permuta 999 veces las etiquetas de especie de la matriz de
rasgos manteniendo fija la matriz de comunidad:

- **SES < 0** — las especies que coexisten son *más similares* que un sorteo
  aleatorio de igual riqueza: agrupamiento funcional, la firma esperada del
  filtrado urbano.
- **SES > 0** — *más distintas* de lo esperado: sobredispersión.
- **SES ≈ 0** — la zona es un subconjunto aleatorio, y la diferencia entre zonas
  es una diferencia de riqueza y nada más.

El pool es el listado **del campus**, así que la prueba no puede responder si el
campus mismo está filtrado respecto de la flora regional del bosque seco. Eso
requiere un pool regional que este estudio no muestreó.

## Hallazgos e implicaciones de manejo

**Conservar prioritariamente la Zona 4 y su arbolado de gran porte** (*Tabebuia
aurea*, *Gliricidia sepium*, *Mangifera indica*, *Ficus benjamina*). El
argumento no se apoya en los índices funcionales —en aves no discriminan entre
zonas— sino en evidencia convergente: mayor riqueza de flora (26 especies),
mayor cobertura arbórea acumulada (101.3 %), el único estrato de dosel que usan
las aves frugívoras registradas, cinco de las seis especies nativas de fauna, y
la única *Ixora coccinea* del campus, único recurso floral sobre el que se
observó actividad de polinizadores.

**Sembrar por los rasgos que faltan, no por las especies bonitas.** La matriz de
rasgos identifica tres vacíos concretos: el 86 % de las especies son
perennifolias y las seis caducifolias son todas árboles introducidos, de modo que
**no existe ningún arbusto caducifolio**; el 26 % son palmas ornamentales, todas
perennifolias y zoócoras, el bloque más redundante del ensamblaje; y solo 8 de 42
especies están en el extremo conservador del espectro de economía foliar.
Nativas del bosque seco que los llenan: *Handroanthus chrysanthus* (caducifolio,
anemócoro, sustituto funcional de *Tabebuia aurea*), *Cordia lutea* (arbustivo
semicaducifolio, floración prolongada), *Bonellia sprucei* y *Pithecellobium
excelsum*.

**Diversificar el recurso floral**, para que la comunidad de polinizadores no sea
una especie introducida visitando un arbusto introducido.

**Repetir este pipeline al concluir la obra**, con la misma metodología y este
estudio como referencia.

## Limitaciones

- **Tres campañas, una ventana estacional** (transición lluviosa–seca). No se
  captura el recambio estacional completo.
- **El esfuerzo de muestreo varía 4.4 veces** entre zonas (2.9 a 12.9 puntos
  ha⁻¹; cobertura efectiva del 20.5 % al 90.9 %). La rarefacción ayuda con
  conteos; no corrige esto del todo.
- **No hay diseño antes–después ni zona control no intervenida.** Los cambios
  futuros no podrán atribuirse inequívocamente a la remodelación: la línea de
  referencia sostiene descripción, no inferencia causal.
- **n pequeño en todo el estudio** — cuatro zonas, 2 a 26 especies cada una. Poco
  margen para estadística inferencial; los modelos nulos son la herramienta
  adecuada a esta escala.
- **El LDMC está sistemáticamente sobreestimado** (sin rehidratación previa a la
  masa fresca) y las submuestras foliares provienen de un solo individuo por
  especie. El ordenamiento relativo en el eje de economía foliar se mantiene
  —que es lo que usa Gower— pero los valores absolutos no son comparables con
  TRY o GLOPNET. Otras cuatro desviaciones de protocolo están documentadas en
  [`data/README.md`](data/README.md).
- **La Zona 2, el 57.3 % del campus, nunca se evaluó** — ya estaba en obra. La
  línea de referencia cubre el 42.7 % del sitio.
- **Puntos de verificación abiertos**, incluida la discrepancia de redundancia sin
  resolver: [`docs/verification.md`](docs/verification.md).

## Cómo citar

Vea [`CITATION.cff`](CITATION.cff) — GitHub lo renderiza como el botón "Cite this
repository".

## Licencia

Código (`R/`, `docs/`): **MIT**. Datos (`data/`, `Tablas_muestreos.ods`):
**CC BY 4.0**. Vea [`LICENSE`](LICENSE).
