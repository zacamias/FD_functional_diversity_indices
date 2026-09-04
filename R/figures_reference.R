#!/usr/bin/env Rscript
# ==============================================================================
#  FIGURES FOR THE README
# ------------------------------------------------------------------------------
#  Renders the figures shown in the repository README into figures/.
#
#  Two of them (the FRic/FEve/FDiv/FDis heat maps) are drawn from the PUBLISHED
#  values of the thesis, Tables 8 and 9, stored in
#  data/reference_results_thesis.csv. They are the reference the pipeline is
#  meant to reproduce, so the README shows the published numbers rather than a
#  re-run that might silently differ.
#
#  The other two (redundancy and the taxonomic-vs-functional contrast) are
#  COMPUTED here from data/, because they need only a Gower distance and can
#  therefore run without the FD package -- useful on machines where FD is not
#  installable.
#
#  Requires: ggplot2, reshape2, cluster, vegan. (Not FD.)
#  Run from the repository root:  Rscript R/figures_reference.R
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(reshape2); library(cluster); library(vegan)
})
source(file.path("R", "redundancy.R"))
source(file.path("R", "taxonomic_diversity.R"))

dir.create("figures", showWarnings = FALSE)
zone_lab <- function(x) gsub("Zona_", "Zone ", x)

theme_fd <- function() {
  theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(colour = "grey35", size = 9),
          plot.caption = element_text(colour = "grey45", size = 7.5, hjust = 0))
}

fmt_index <- function(value, index) {
  ifelse(is.na(value), "NA",
         ifelse(index == "FRic", formatC(value, format = "e", digits = 2),
                formatC(value, format = "f", digits = 3)))
}

## ---- 1. Heat maps from the published thesis values ---------------------------
ref <- read.csv("data/reference_results_thesis.csv", stringsAsFactors = FALSE)

heat_group <- function(g, title, subtitle, outfile) {
  d <- ref[ref$grupo == g, c("Zona", "FRic", "FEve", "FDiv", "FDis")]
  d$Zona <- zone_lab(d$Zona)
  # z-score per index column; an all-NA or constant column stays at 0.
  z <- d
  for (cc in c("FRic", "FEve", "FDiv", "FDis")) {
    v <- d[[cc]]; s <- sd(v, na.rm = TRUE)
    z[[cc]] <- if (!is.na(s) && s > 0) (v - mean(v, na.rm = TRUE)) / s else rep(0, length(v))
  }
  lz <- melt(z, id.vars = "Zona", variable.name = "Index", value.name = "z")
  lv <- melt(d, id.vars = "Zona", variable.name = "Index", value.name = "value")
  la <- merge(lz, lv, by = c("Zona", "Index"))
  la$lab <- fmt_index(la$value, la$Index)
  la$z[is.na(la$value)] <- NA

  g1 <- ggplot(la, aes(x = Index, y = Zona, fill = z)) +
    geom_tile(colour = "white", linewidth = 0.7) +
    geom_text(aes(label = lab), size = 3.2) +
    scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                         midpoint = 0, na.value = "grey88", name = "z-score") +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL,
         caption = "Colour is standardized WITHIN each column, so it compares zones for one index only.\nIt is not comparable across columns. Grey = index not computable (see index_feasibility_*.csv).") +
    theme_fd()
  ggsave(file.path("figures", paste0(outfile, ".png")), g1, width = 7.2, height = 3.9, dpi = 200)
  invisible(g1)
}

heat_group("flora",
           "Functional diversity of flora, by zone",
           "Las Penas campus (ESPOL). FRic is a volume (note the scientific notation); FEve, FDiv and FDis are 0-1.\nZone 1 holds only 2 plant species, below the 3 needed for a convex hull.",
           "heatmap_flora")

heat_group("aves",
           "Functional diversity of birds, by zone",
           "Campaigns pooled. FRic is NA because zones 3, 4 and 5 share exactly the same four species,\nso an unweighted volume cannot discriminate between them.",
           "heatmap_birds")

## ---- 2. Functional redundancy (computed from data/) --------------------------
gower_dist <- function(traits) {
  # cluster::daisy with metric="gower" reproduces FD::gowdis for mixed traits
  # (both are Gower 1971 with range rescaling of the numeric variables).
  as.matrix(cluster::daisy(traits, metric = "gower"))
}

prep_traits <- function(df, continuous, categorical, log_cols = character(0)) {
  rownames(df) <- df$Nombre_cientifico
  df$Nombre_cientifico <- NULL
  for (cc in log_cols)    df[[cc]] <- log10(df[[cc]])
  for (ct in categorical) df[[ct]] <- as.factor(df[[ct]])
  df[, c(continuous, categorical), drop = FALSE]
}

prep_comm <- function(abund, species_ref) {
  rownames(abund) <- abund$Nombre_cientifico
  abund$Nombre_cientifico <- NULL
  m <- t(as.matrix(abund))
  m <- m[, colnames(m) %in% species_ref, drop = FALSE]
  m <- m[rowSums(m) > 0, , drop = FALSE]
  m[, colSums(m) > 0, drop = FALSE]
}

rd <- function(f) read.csv(file.path("data", f), stringsAsFactors = FALSE)

# --- flora
tf <- prep_traits(rd("traits_flora.csv"),
                  continuous  = c("H_m", "LDMC", "SLA"),
                  categorical = c("forma_crecimiento","fenologia","forma_hoja","dispersion"),
                  log_cols    = c("H_m", "SLA"))
cf <- prep_comm(rd("abundance_flora.csv"), rownames(tf))
tf <- tf[colnames(cf), , drop = FALSE]
red_flora <- functional_redundancy(cf, gower_dist(tf)); red_flora$grupo <- "Flora"
tax_flora <- taxonomic_diversity(cf, counts = FALSE)

# --- birds, campaigns pooled
meta  <- rd("species_metadata.csv")
birds <- meta$Nombre_cientifico[meta$grupo == "Ave"]
ab <- lapply(1:3, function(i) {
  a <- rd(sprintf("abundance_fauna_survey%d.csv", i))
  a[a$Nombre_cientifico %in% birds, , drop = FALSE]
})
zones <- setdiff(names(ab[[1]]), "Nombre_cientifico")
all_sp <- sort(unique(unlist(lapply(ab, `[[`, "Nombre_cientifico"))))
ab_birds <- data.frame(Nombre_cientifico = all_sp)
for (z in zones) {
  ab_birds[[z]] <- sapply(all_sp, function(s)
    sum(sapply(ab, function(a) { v <- a[a$Nombre_cientifico == s, z]; if (length(v)) v else 0 })))
}
tb_raw <- do.call(rbind, lapply(1:3, function(i) rd(sprintf("traits_fauna_survey%d.csv", i))))
tb_raw <- tb_raw[!duplicated(tb_raw$Nombre_cientifico) & tb_raw$Nombre_cientifico %in% birds, ]
tb <- prep_traits(tb_raw, continuous = "masa_g",
                  categorical = c("gremio_trofico","estrato"), log_cols = "masa_g")
cb <- prep_comm(ab_birds, rownames(tb))
tb <- tb[colnames(cb), , drop = FALSE]
red_birds <- functional_redundancy(cb, gower_dist(tb)); red_birds$grupo <- "Birds"
tax_birds <- taxonomic_diversity(cb, counts = TRUE)

red <- rbind(red_flora, red_birds)
red$Zona <- zone_lab(red$Zone)
write.csv(red, file.path("figures", "redundancy_values.csv"), row.names = FALSE)

# The thesis reports relative redundancy for birds in Table 9 and, for flora,
# only as a range in the text (0.867-0.883). Plotting the published values next
# to this repository's recomputation makes the comparison explicit instead of
# silently replacing one with the other. THEY CURRENTLY DISAGREE -- see
# docs/verification.md, open item V1.
ref_fr <- ref[!is.na(ref$FR_rel), c("grupo", "Zona", "FR_rel")]
ref_fr$grupo <- ifelse(ref_fr$grupo == "aves", "Birds", "Flora")
ref_fr$Zona  <- zone_lab(ref_fr$Zona)
ref_fr$fuente <- "Thesis (Table 9)"

rep_fr <- red[!is.na(red$FR_rel), c("grupo", "Zona", "FR_rel")]
rep_fr$fuente <- "This repository (Gower)"

fr_all <- rbind(rep_fr, ref_fr)
fr_all$fuente <- factor(fr_all$fuente,
                        levels = c("This repository (Gower)", "Thesis (Table 9)"))

g_red <- ggplot(fr_all, aes(x = Zona, y = FR_rel, fill = fuente)) +
  geom_col(position = position_dodge(width = 0.78, preserve = "single"), width = 0.68) +
  geom_text(aes(label = formatC(FR_rel, format = "f", digits = 3)),
            position = position_dodge(width = 0.78, preserve = "single"),
            vjust = -0.4, size = 2.9) +
  facet_wrap(~ grupo) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.14))) +
  scale_fill_manual(values = c("This repository (Gower)" = "#2E7D5B",
                               "Thesis (Table 9)" = "#C08A2E"), name = NULL) +
  labs(title = "Functional redundancy by zone",
       subtitle = "Ricotta et al. (2016): the share of species diversity that does NOT translate into functional diversity.\nThe two series use the same formula (1 - Q/D) but do not agree - an open verification item, not a settled result.",
       x = NULL, y = "Relative redundancy  (1 - Q/D)",
       caption = "Repository values computed from data/ with a Gower distance. The thesis reports per-zone values only for birds (Table 9);\nfor flora it reports the range 0.867-0.883. Zones with a single species are omitted: redundancy is undefined there, not zero.\nSee docs/verification.md, item V1.") +
  theme_fd() + theme(legend.position = "top")
ggsave("figures/redundancy_by_zone.png", g_red, width = 8, height = 4.4, dpi = 200)

## ---- 3. Taxonomic vs functional diversity ------------------------------------
tax_flora$grupo <- "Flora"; tax_birds$grupo <- "Birds"
tax <- rbind(tax_flora, tax_birds)
tax$Zona <- zone_lab(tax$Zone)
write.csv(tax, file.path("figures", "taxonomic_values.csv"), row.names = FALSE)

tvf <- merge(tax[, c("Zona","grupo","S","Simpson")],
             red[, c("Zona","grupo","RaoQ","FR_rel")],
             by = c("Zona","grupo"))
lg <- melt(tvf[, c("Zona","grupo","Simpson","RaoQ")],
           id.vars = c("Zona","grupo"), variable.name = "Measure", value.name = "value")
lg$Measure <- factor(lg$Measure, levels = c("Simpson","RaoQ"),
                     labels = c("Species diversity (Simpson D)",
                                "Functional diversity (Rao's Q)"))

g_tvf <- ggplot(lg, aes(x = Zona, y = value, fill = Measure)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  facet_wrap(~ grupo) +
  scale_fill_manual(values = c("#8C9EA8", "#08519C"), name = NULL) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.08))) +
  labs(title = "Species diversity is not functional diversity",
       subtitle = "The gap between the two bars IS the functional redundancy: species present whose role is already covered.",
       x = NULL, y = "Index value",
       caption = "Simpson D and Rao's Q share the same 0-1 scale and are directly comparable; Q uses the Gower trait distance.") +
  theme_fd() + theme(legend.position = "top")
ggsave("figures/taxonomic_vs_functional.png", g_tvf, width = 8, height = 4.2, dpi = 200)

cat("\nFigures written to figures/\n")
print(red[, c("grupo","Zona","n_species","Simpson_D","RaoQ","FR_rel")], row.names = FALSE)
