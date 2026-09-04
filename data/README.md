# Data dictionary · Diccionario de variables

All tables in this folder are UTF-8 CSV with a header row. Column names are kept
in Spanish to match the source data dictionary of the field notebooks and of
`Tablas_muestreos.ods`; this file gives the English gloss, the unit and the
provenance of every variable.

Todas las tablas de esta carpeta son CSV en UTF-8 con fila de encabezado. Los
nombres de columna se mantienen en español para coincidir con el diccionario de
datos original de las fichas de campo y de `Tablas_muestreos.ods`.

**Source of truth.** These CSVs were derived from `Tablas_muestreos.ods`
(9 sheets), cleaned and harmonized. The `.ods` is kept in the repository as the
raw field record; the CSVs are the analysis-ready form and are what the pipeline
reads.

---

## 1. Study design

Las Peñas campus (ESPOL, Guayaquil, Ecuador; ~2°11′S, 79°53′W), 3.25 ha.
Five zones; **Zone 2 was excluded** because it was under civil works during the
study. The four assessed zones total 13 893 m² (42.7 % of the campus).

- **Flora:** a single complete census (no sampling — every vascular plant
  individual in each zone was recorded).
- **Fauna:** three monthly campaigns (May–July 2026), covering the transition
  from the wet to the dry season. Birds by fixed-radius (15 m) point counts,
  reptiles and pollinators by timed 15-minute walks.

### `zones.csv`

| Column | Unit | Description |
|---|---|---|
| `Zona` | — | Zone identifier. Matches the column names of the abundance matrices. |
| `area_m2` | m² | Total zone area. |
| `pct_campus` | % | Share of the 3.25 ha campus. |
| `puntos_conteo` | count | Number of bird point counts placed in the zone. |
| `puntos_ha` | points ha⁻¹ | Sampling effort per unit area. **Varies 4.4× between zones** (2.9 in Zone 3 to 12.9 in Zone 5) — this must be considered when comparing zones. |
| `superficie_cubierta_pct` | % | Share of the zone effectively covered by the point counts (706.9 m² per point, 15 m radius). |
| `detecciones_fauna` | count | Total detections pooled over the three campaigns and all three faunal groups. |
| `especies_fauna` | count | Cumulative species richness of fauna. |
| `evaluada` | logical | Whether the zone was assessed. |

---

## 2. Flora

### `traits_flora.csv` — species × trait matrix (42 species)

| Column | Unit | Type | Description |
|---|---|---|---|
| `Nombre_cientifico` | — | key | Scientific name, verified in Plants of the World Online (POWO, 2026), cross-checked in World Flora Online. Joins to `abundance_flora.csv` and `species_metadata.csv`. |
| `forma_crecimiento` | — | nominal | Growth form: `Árbol`, `Árbol pequeño`, `Arbusto`, `Arbusto suculento`, `Cactus`, `Cícada`, `Hierba`, `Palma`, `Suculenta`, `Trepadora`. From specialized botanical/horticultural literature and species datasheets. |
| `H_m` | **m** | continuous | Maximum height of the tallest individual of the species in the campus. Small individuals measured with tape; trees with a percentage clinometer, applying `H = d × (%crown − %base) / 100`. `Spathodea campanulata` (22.35 m) comes from a direct field measurement (the cell was blank in the original `.ods`). |
| `fenologia` | — | nominal | Leaf phenology: `Perenne` (evergreen) or `Caducifolia` (deciduous). |
| `forma_hoja` | — | nominal | Leaf form: `Simple`, `Compuesta`, `Compuesta (palma)`, `Suculenta/reducida`. |
| `dispersion` | — | nominal | Dispersal syndrome: `Zoocoria`, `Anemocoria`, `Autocoria`, `Hidrocoria`. |
| `LDMC` | **g g⁻¹** (dimensionless) | continuous | Leaf dry matter content = dry mass / fresh mass. **Note the unit:** the thesis reports LDMC in mg g⁻¹; these values are the same quantity expressed as a fraction. Multiply by 1000 to recover the thesis figures (e.g. `0.338` here = 338.2 mg g⁻¹ in Table 4). |
| `SLA` | **mm² g⁻¹** | continuous | Specific leaf area = leaf area / dry mass. **Note the unit:** the thesis reports SLA in mm² mg⁻¹. Divide by 1000 to recover the thesis figures (e.g. `6329.14` here = 6.33 mm² mg⁻¹ in Table 4). |

**Laboratory protocol for LDMC and SLA** (Pérez-Harguindeguy et al., 2013):
three leaf samples per species (42 species, 126 samples) from the largest
individual; fresh mass on an analytical balance (0.0001 g); leaf area digitized
by calibrated digital photography and measured in ImageJ; oven-dried at 70 °C
for 48 h, cooled in a silica-gel desiccator, then dry mass recorded.

> **Known deviations from the reference protocol** (documented in §2.6 of the
> thesis, and relevant to anyone reusing these values):
> 1. **No rehydration** before fresh-mass determination, and collection took
>    place 11:00–13:00, at peak evaporative demand. LDMC is therefore
>    **systematically overestimated**. Because the bias acts in the same
>    direction for every species, the relative ordering along the leaf-economics
>    axis — which is what the Gower distance uses — remains interpretable, but
>    the absolute values are **not directly comparable to international trait
>    databases** (TRY, GLOPNET).
> 2. The three leaves per species came from **one individual**, so they are
>    subsamples, not independent biological replicates. Several ornamental
>    species were represented by one or two individuals campus-wide.
> 3. The largest individual was chosen, biasing toward **sun leaves** (lower
>    SLA). Applied consistently across species.
> 4. Leaf area from photography rather than flatbed scanning; perspective
>    distortion mitigated by a parallel camera plane and an in-frame scale.
> 5. Drying limited to **48 h instead of 72 h** (shared oven). Residual moisture
>    cannot be excluded in the succulent and aphyllous species, which would
>    overestimate SLA and underestimate LDMC **for that group only**.

### `abundance_flora.csv` — species × zone abundance

| Column | Unit | Description |
|---|---|---|
| `Nombre_cientifico` | — | Joins to `traits_flora.csv`. |
| `Zona_4`, `Zona_5`, `Zona_3`, `Zona_1` | **fraction 0–1** | Cover: horizontal crown projection of all individuals of the species divided by the total zone area. **Note the unit:** the thesis reports this as a percentage. Multiply by 100 (e.g. `0.6132` = 61.32 %). `0` means the species is absent from the zone. |

> Cumulative cover exceeds 100 % in Zone 4 (101.3 %) and Zone 5 (110.0 %)
> because crowns of different strata overlap vertically — expected when
> quantifying crown projection rather than ground cover.

---

## 3. Fauna

Three campaigns are stored as three separate pairs of files
(`traits_fauna_survey{1,2,3}.csv` / `abundance_fauna_survey{1,2,3}.csv`) because
the species pool grew between campaigns: 6 species in campaign 1, 9 in
campaign 2, 8 in campaign 3. Each trait file contains exactly the species
detected in that campaign.

### `traits_fauna_survey{1,2,3}.csv`

| Column | Unit | Type | Description |
|---|---|---|---|
| `Nombre_cientifico` | — | key | Nomenclature: eBird/Clements (Clements et al., 2025) for birds; Reptiles del Ecuador (Torres-Carvajal et al., 2026) for reptiles; iNaturalist used for photographic contrast only, not as a nomenclatural authority. |
| `masa_g` | **g** | continuous | Body mass. Birds from EltonTraits 1.0 (Wilman et al., 2014) and Ridgely & Greenfield (2001); reptiles from Cruz-García et al. (2026) and herpetological literature; *Apis mellifera* from entomological literature. Where the source gives a range, the **midpoint** is used (e.g. *Columbina buckleyi* 65–72 g → 68.5 g). |
| `gremio_trofico` | — | nominal | Trophic guild: `Granívoro`, `Frugívoro`, `Omnívoro`, `Insectívoro`, `Herbívoro`, `Nectarívoro`. |
| `estrato` | — | nominal | Activity stratum: `Terrestre`, `Sotobosque`, `Dosel`, `Aéreo`. |

Body-mass values actually used:

| Species | `masa_g` | Source range | Note |
|---|---|---|---|
| *Apis mellifera* | 0.1 | 0.09–0.12 g | |
| *Anolis sagrei* | 6.0 | 4–8 g | Males ≈8 g, females ≈4 g |
| *Troglodytes musculus* | 12.0 | 12 g | |
| *Euphonia xanthogaster* | 12.5 | 9–16 g | |
| *Thraupis episcopus* | 36.0 | 27–45 g | |
| *Columbina buckleyi* | 68.5 | 65–72 g | |
| *Quiscalus mexicanus* | 155.0 | 155 g | |
| *Columba livia* | 370.0 | 370 g | |
| *Iguana iguana* | 2000.0 | 1 500–10 000 g | Lower end used: campus individuals are subadult |

### `abundance_fauna_survey{1,2,3}.csv`

| Column | Unit | Description |
|---|---|---|
| `Nombre_cientifico` | — | Joins to the trait file of the same campaign. |
| `Zona_4`, `Zona_5`, `Zona_3`, `Zona_1` | **count of detections** | Number of detections of the species in the zone during that campaign. |

> **These are detections, not abundances.** For birds and reptiles,
> independence between successive observations cannot be guaranteed, so the
> values are reported as detections per species and zone.
>
> **The *Apis mellifera* values are not direct counts.** They are extrapolated
> estimates from four 2.5 m linear sampling units placed in four different
> *Ixora coccinea* planters (10 m of floral resource per campaign), scaled to
> the 304 linear metres the species occupies campus-wide (extrapolation factor
> 30.4). The 95 % confidence intervals are wide: campaign 1, 243 (85–401);
> campaign 2, 365 (207–523); campaign 3, 365 (91–638). The units were placed in
> the four **longest** planters and not chosen at random; because *A. mellifera*
> concentrates foragers on the richest patches through dance recruitment
> (Beekman & Ratnieks, 2000), the extrapolated figure should be read as an
> **upper bound**, and as an index of visitation intensity per unit of floral
> resource — **not** as a population size estimate. Consequently it is **not on
> the same scale as the vertebrate detections**, which is why the pipeline also
> runs a birds-only analysis (see `R/functional_diversity_analysis.R`).

---

## 4. `species_metadata.csv`

| Column | Description |
|---|---|
| `Nombre_cientifico` | Joins to every other table. |
| `grupo` | `Flora`, `Ave`, `Reptil`, `Insecto`. |
| `origen` | `Nativa` or `Introducida`. All 42 plant species are introduced ornamentals; of the 9 animal species, 6 are native. |
| `notas` | Free text: why the record matters, or how a value was obtained. |

Not used in the distance computation — this table is provenance and
interpretation metadata, deliberately kept out of the trait matrix so it cannot
leak into the Gower distance.

---

## 5. Zone column order

The abundance matrices store zones in the order `Zona_4, Zona_5, Zona_3,
Zona_1`, inherited from the original spreadsheet. The pipeline transposes and
reorders them, so this order carries no meaning; do not read it as a ranking.

---

## 6. Missing data

There are no `NA` cells in these tables. A `0` in an abundance matrix means a
true absence of the species from the zone, not a missing observation. `NA`
appears only in the **results**, where it flags an index that is not computable
for a given zone (see the richness diagnostic in the pipeline output).

---

## 7. Reproducing the thesis tables from these files

| Thesis table | File | Transformation |
|---|---|---|
| Table 4 (flora traits) | `traits_flora.csv` | `LDMC × 1000` → mg g⁻¹; `SLA / 1000` → mm² mg⁻¹ |
| Table 5 (flora cover) | `abundance_flora.csv` | `× 100` → % |
| Table 6 (fauna traits) | `traits_fauna_survey2.csv` + `species_metadata.csv` | Campaign 2 holds the full 9-species pool |
| Table 7 (fauna detections) | `abundance_fauna_survey{1,2,3}.csv` | Bind the three campaigns |
| Tables 8–9 (indices) | generated | `Rscript R/functional_diversity_analysis.R` |
| Table 10 (reptiles/pollinators) | `abundance_fauna_survey{1,2,3}.csv` | Filter `grupo` ∈ {Reptil, Insecto} via `species_metadata.csv` |

---

## References for the trait sources

- Beekman, M., & Ratnieks, F. L. W. (2000). Long-range foraging by the honey-bee, *Apis mellifera* L. *Functional Ecology*, 14(4), 490–496.
- Clements, J. F., et al. (2025). *The eBird/Clements checklist of birds of the world*.
- Cruz-García, et al. (2026). Urban herpetofauna of Guayaquil.
- Pérez-Harguindeguy, N., et al. (2013). New handbook for standardised measurement of plant functional traits worldwide. *Australian Journal of Botany*, 61(3), 167–234.
- POWO (2026). *Plants of the World Online*. Royal Botanic Gardens, Kew.
- Ridgely, R. S., & Greenfield, P. J. (2001). *The Birds of Ecuador*. Cornell University Press.
- Taylor, L. R. (1961). Aggregation, variance and the mean. *Nature*, 189, 732–735.
- Torres-Carvajal, O., et al. (2026). *Reptiles del Ecuador*.
- Wilman, H., et al. (2014). EltonTraits 1.0. *Ecology*, 95(7), 2027.
