# Verification log

This file records every point where the code, the data and the thesis
(*Caracterización de la diversidad funcional del campus Las Peñas durante su
remodelación*, ESPOL, 2026) were compared, what was found, and what was done
about it.

It exists because a research compendium is only worth something if the reader
can tell **which numbers have been checked** and which have not. Items marked
**RESOLVED** were fixed in this repository; items marked **OPEN** are real
discrepancies that still need the author's decision, and they are stated here
rather than papered over.

The thesis is treated as the authoritative record throughout: where the code and
the thesis disagreed on a fact, the code was changed, not the thesis.

---

## OPEN items

### V1 — Functional redundancy values do not reproduce

**Status: OPEN — this is the most important item in this file.**

The thesis reports relative functional redundancy (Ricotta et al. 2016) of
**0.794–0.872** for birds (Table 9) and **0.867–0.883** for flora (text, §3.2),
and builds the central resilience argument on them: "about 87 % of species
diversity does not translate into functional diversity."

An independent implementation of the same formula in `R/redundancy.R`, applied
to the same data, gives:

| Group | Zone | This repository | Thesis |
|---|---|---|---|
| Birds | Zone 1 | 0.406 | 0.794 |
| Birds | Zone 3 | 0.206 | 0.829 |
| Birds | Zone 4 | 0.304 | 0.872 |
| Birds | Zone 5 | 0.253 | 0.850 |
| Flora | Zone 1 | 0.540 | 0.867–0.883 (range) |
| Flora | Zone 3 | 0.530 | " |
| Flora | Zone 4 | 0.495 | " |
| Flora | Zone 5 | 0.574 | " |

Both use `FR_rel = 1 − Q/D`, with `D` the Gini–Simpson index and `Q` Rao's
quadratic entropy. The gap is large and systematic (the thesis values are
roughly 1.5–4× higher), which points at the **dissimilarity used inside Q**
rather than at an arithmetic slip. The likely candidates, in order:

1. **Distance scaling.** Ricotta et al. require `d ∈ [0, 1]`. If the thesis
   computed `Q` on a distance already rescaled by something other than the Gower
   maximum — for example on PCoA axis coordinates, or on the Cailliez-corrected
   distance divided by its own maximum — `Q` shrinks and `FR_rel` rises.
2. **Squared vs unsquared distances.** Using `d²` instead of `d` in
   `Q = Σ Σ d_ij p_i p_j` lowers `Q` substantially for `d < 1`, which would push
   the values in exactly the observed direction.
3. **Which distance.** Gower on the raw traits (what this repository uses) vs
   Gower after the `log10` transformations vs a Euclidean distance in the
   corrected PCoA space.

**What is needed:** the lines of the original script that computed `D` and `Q`.
Once the definition is settled, either `R/redundancy.R` is corrected to match
it, or the thesis figure is revised. Until then, `figures/redundancy_by_zone.png`
plots both series side by side and labels them, so the reader is never shown one
without the other.

This matters beyond bookkeeping: the "87 % redundancy" claim is what supports
the buffering-capacity conclusion, so it should be the best-documented number in
the study.

### V2 — Internal inconsistencies within the thesis

Section 3.2 and Tables 8–9 disagree with each other in four places. The tables
are treated as primary here and are what `data/reference_results_thesis.csv`
stores, but the text should be reconciled:

| Quantity | Table 8/9 | Text §3.2 |
|---|---|---|
| Flora FEve, Zone 3 | 0.560 | 0.648 |
| Flora FEve, Zone 5 | 0.371 | 0.366 |
| Flora FDis, Zone 3 | 0.069 | 0.065 |
| Flora FDis, Zone 5 | 0.293 | 0.279 |

Separately, §3.2 states that Zone 1 birds reach "a functional richness 2.5 times
higher" than the other zones, while Table 9 reports **FRic = NA for every bird
zone**. Both statements cannot stand. The most likely reading is that a FRic
value was computed for Zone 1 during exploration and then dropped from the table
when FRic was judged uninformative; if so, the sentence needs rewriting, because
as it stands it cites a number the tables say does not exist.

### V3 — Trait category schemes differ between Table 4 and the working data

`data/traits_flora.csv` and Table 4 of the thesis classify the same species
differently in two columns:

| Species | Column | Working data | Thesis Table 4 |
|---|---|---|---|
| *Beaucarnea recurvata* | growth form / leaf form | `Árbol` / `Simple` | `Suculenta/Cactus` / `Lineal` |
| *Cycas revoluta* | growth form | `Cícada` | `Palma` |
| *Morinda citrifolia* | growth form | `Árbol pequeño` | `Arbusto` |
| *Cordia sebestena* | phenology | `Perenne` | `Caducifolia` |
| *Manihot esculenta* | leaf form | `Simple` | `Palmada` |
| *Tabebuia aurea* | leaf form | `Compuesta` | `Palmada` |

The working data uses a finer growth-form scheme (10 levels, splitting cacti,
succulents and cycads) and a coarser leaf-form scheme (4 levels, collapsing
`Palmada` and `Lineal`). **These are not typos, they are two different
classifications**, and because Gower weights every categorical trait equally,
the choice changes the distance matrix and therefore every index.

The working-data scheme was kept, because it is internally consistent and is
what the pipeline was written against. But *Cordia sebestena*'s phenology
(evergreen vs deciduous) is a factual disagreement, not a scheme difference, and
should be resolved from the field notes.

Note also that Table 4 as typeset carries eight column headers over seven
columns of data: the height column appears to have been lost in layout, so the
values printed under "Altura (m)" are LDMC and those under "LDMC (mg/g)" are
SLA. The heights in `data/traits_flora.csv` are the field measurements and are
not affected.

---

## RESOLVED items

### V4 — Species misidentification in the code: *Claravis pretiosa* → *Columbina buckleyi*

The analysis script named the ground dove **`Claravis pretiosa`** (Blue Ground
Dove). Every table in the thesis names **`Columbina buckleyi`** (Ecuadorian
Ground Dove). The body mass in the script, 68.5 g, is the midpoint of the 65–72 g
range the thesis gives for *C. buckleyi*; *C. pretiosa* is a smaller, forest-
interior species that would be a remarkable record on an urban campus.

This was a name error in the code, not a different bird. **Corrected** in all
four data files.

### V5 — Body masses realigned to Table 6

| Species | Was | Now | Source range (Table 6) |
|---|---|---|---|
| *Columba livia* | 340.0 g | **370.0 g** | 370 g |
| *Quiscalus mexicanus* | 180.0 g | **155.0 g** | 155 g |
| *Euphonia xanthogaster* | 13.0 g | **12.5 g** | 9–16 g (midpoint) |

*Iguana iguana* was left at 2000 g. Table 6 gives 1 500–10 000 g, whose midpoint
(5 750 g) would be a large adult; campus individuals are subadult, so the lower
end is the defensible choice. It is now stated explicitly in `data/README.md`
instead of being an unexplained constant.

### V6 — *Apis mellifera*, campaign 3: 360 → 365 detections

The abundance matrix for campaign 3 carried 360; Tables 7 and 10 both give 365.
**Corrected.**

### V7 — PCoA correction: Lingoes → Cailliez

The script defaulted to `lingoes`; the thesis states Cailliez (Table 8 note) and
so did the repository description. Lingoes adds a smaller constant and distorts
the distance structure less, so it is arguably the better choice on the merits —
but the published tables are Cailliez values, and a repository whose purpose is
to reproduce them must default to what produced them.

**Default changed to `cailliez`**, with the trade-off documented in the script
header. `pcoa_correction <- "lingoes"` remains available as a one-line
sensitivity check, which is how the comparison should be reported.

### V8 — No richness control

Every index here scales with richness, and the assessed zones span 2 to 26 plant
species. Comparing zones directly therefore confounds "how many species" with
"how they are arranged in trait space".

**Added** `R/null_models.R`: 999 trait-label permutations produce a standardized
effect size and an exact two-tailed p-value per zone and index. The question it
answers is the one that was previously unanswerable — *is this zone's assemblage
more functionally clustered than an equally rich random draw from the campus
pool?* Note the pool is campus-level, so the result cannot speak to whether the
campus as a whole is filtered relative to the regional dry-forest flora.

### V9 — Redundancy not computed by the pipeline

The thesis reports redundancy; the script did not compute it. **Added**
`R/redundancy.R` (see V1 for the open disagreement about its value).

### V10 — Undocumented NA values

FRic/FEve/FDiv are NA for Zone 1 flora with no machine-readable explanation.
**Added** a feasibility diagnostic (`index_feasibility_*.csv`) that states, per
zone, which indices are computable and why not, so an NA is a documented
structural limit rather than a gap.

### V11 — Unit ambiguity in LDMC and SLA

The data carry LDMC as a dimensionless fraction (g g⁻¹) and SLA in mm² g⁻¹,
while the thesis reports mg g⁻¹ and mm² mg⁻¹. The values are the same quantity,
a factor of 1000 apart, but nothing said so. **Documented** in
`data/README.md`, with the conversion needed to recover each thesis table.

---

## How to re-run the checks

```bash
Rscript R/figures_reference.R          # redundancy + taxonomic recomputation (no FD needed)
Rscript R/functional_diversity_analysis.R   # full pipeline (needs FD)
```

Then compare `results/*/indices_dbFD_*.csv` against
`data/reference_results_thesis.csv`. Any zone where the two differ by more than
rounding is a new item for this file.
