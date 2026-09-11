#!/usr/bin/env Rscript
# ==============================================================================
#  FUNCTIONAL DIVERSITY ANALYSIS — Las Penas Campus (ESPOL)
#  Capstone Project · Miguel Angel Moran Piedrahita · 2026
# ------------------------------------------------------------------------------
#  Goal: compute FRic, FEve, FDiv, FDis, Rao's Q and functional redundancy per
#        ZONE for flora and for fauna, and produce heat maps.
#
#  Methodological decisions:
#   - Flora abundance = crown-projection cover (fraction 0-1); fauna = detection
#     counts. See data/README.md for units and their caveats.
#   - Distance        = Gower (mixed traits; internal range rescaling).
#   - log10 on skewed continuous traits: mass (fauna), height and SLA (flora).
#   - Ordination      = PCoA with CAILLIEZ correction for negative eigenvalues
#                       (this is what the thesis used; see note in section 0).
#   - Engines: FD::dbFD, with mFD as an optional cross-check.
#   - Redundancy      = Ricotta et al. (2016), FR = D - Q, reported as 1 - Q/D.
#   - Null models     = 999 trait-label permutations -> richness-independent SES.
#   - Taxonomic diversity computed in parallel, so "few functions despite many
#     species" is shown rather than asserted.
#   - Sensitivity analysis: indices recomputed with raw abundance, log(x+1) and
#     presence/absence, to separate the "dominance" effect (Apis mellifera,
#     90.1-95.7 % of Zone 4 fauna records) from the "composition" effect.
#   - BIRDS-ONLY analysis pooled across the three campaigns, which is the
#     analysis the thesis reports (Table 9), run alongside the all-fauna
#     per-campaign analysis that motivates it.
#
#  Data live in data/ as CSV. See data/README.md for the full data dictionary,
#  units, laboratory protocol deviations, and how to reproduce each thesis table.
# ==============================================================================

## ---- 0. CONFIGURATION --------------------------------------------------------
dir_data        <- "data"
dir_base        <- "results"      # all outputs are written here
pcoa_correction <- "cailliez"     # see note below
use_log_height  <- TRUE           # log10 on flora height (skew ~1.16)
n_perm          <- 999            # permutations for the null models
run_null_models <- TRUE           # set FALSE for a fast run (this is the slow part)
set.seed(123)

#  NOTE ON THE PCoA CORRECTION.
#  Gower distances on mixed traits are generally NOT Euclidean, which produces
#  negative eigenvalues in the PCoA and makes the convex hull behind FRic
#  ill-defined. Two standard corrections exist:
#    - Cailliez (1983): adds a constant to the DISTANCES.
#    - Lingoes (1971):  adds a constant to the SQUARED distances; the constant
#                       is smaller, so the distance structure is distorted less.
#  Lingoes is the less invasive of the two and is often preferable on those
#  grounds. Cailliez is nevertheless the default here for one reason: it is what
#  the thesis used, and the published Tables 8-9 are Cailliez values. Keeping
#  the default aligned with the published record means this repository
#  reproduces the thesis rather than quietly disagreeing with it.
#  Set pcoa_correction <- "lingoes" to check how much the choice actually
#  matters; do so as a documented sensitivity check, not as a silent change.

## ---- 1. PACKAGES -------------------------------------------------------------
pkgs_required <- c("FD", "ape", "ade4", "ggplot2", "reshape2", "vegan")
missing <- pkgs_required[
  !vapply(pkgs_required, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing) > 0) {
  stop(
    "Missing required packages: ", paste(missing, collapse = ", "), ".\n",
    "This project pins its environment with renv. Restore it with:\n",
    "  install.packages('renv'); renv::restore()\n",
    "Installing these manually will not reproduce the published results.",
    call. = FALSE
  )
}
invisible(lapply(pkgs_required, function(p) {
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}))
# OPTIONAL: mFD. If missing or not installable (it needs cmake, libgmp-dev,
# libnlopt-cxx-dev), the script continues with FD::dbFD only and says so.
mFD_available <- requireNamespace("mFD", quietly = TRUE)
if (mFD_available) suppressPackageStartupMessages(library(mFD))

# Project modules
source(file.path("R", "redundancy.R"))
source(file.path("R", "null_models.R"))
source(file.path("R", "taxonomic_diversity.R"))

## ---- 2. DIRECTORIES ----------------------------------------------------------
dir_flora <- file.path(dir_base, "flora")
dirs_fau  <- file.path(dir_base, "fauna", c("survey_1", "survey_2", "survey_3"))
dir_birds <- file.path(dir_base, "birds_pooled")
dir_fcomp <- file.path(dir_base, "fauna", "cross_survey_comparison")
dir_logs  <- file.path(dir_base, "_logs")
for (d in c(dir_flora, dirs_fau, dir_birds, dir_fcomp, dir_logs)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
log_file <- file.path(dir_logs, paste0("log_", format(Sys.Date(), "%Y%m%d"), ".txt"))
log_msg  <- function(...) {
  line <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ...)
  cat(line, "\n"); cat(line, "\n", file = log_file, append = TRUE)
}
log_msg("Analysis start. R ", as.character(getRversion()),
        " | mFD available: ", mFD_available,
        " | PCoA correction: ", pcoa_correction)

## ---- 3. DATA -----------------------------------------------------------------
read_data <- function(f) {
  path <- file.path(dir_data, f)
  if (!file.exists(path))
    stop("Missing data file: ", path,
         "\nRun this script from the repository root, not from R/.")
  read.csv(path, stringsAsFactors = FALSE, encoding = "UTF-8",
           check.names = FALSE)
}

traits_flora <- read_data("traits_flora.csv")
abund_flora  <- read_data("abundance_flora.csv")
traits_fau   <- lapply(1:3, function(i) read_data(sprintf("traits_fauna_survey%d.csv", i)))
abund_fau    <- lapply(1:3, function(i) read_data(sprintf("abundance_fauna_survey%d.csv", i)))
sp_meta      <- read_data("species_metadata.csv")

## ---- 4. HELPER FUNCTIONS -----------------------------------------------------

# 4.1 TRAIT matrix: rownames = species, log10 on skewed continuous, factors.
prep_traits <- function(df, continuous, categorical, log_cols = character(0)) {
  rownames(df) <- df$Nombre_cientifico
  df$Nombre_cientifico <- NULL
  for (cc in log_cols)    df[[cc]] <- log10(df[[cc]])
  for (ct in categorical) df[[ct]] <- as.factor(df[[ct]])
  df[, c(continuous, categorical), drop = FALSE]
}

# 4.2 COMMUNITY matrix: zones (rows) x species (columns).
prep_comm <- function(abund, species_ref, group) {
  rownames(abund) <- abund$Nombre_cientifico
  abund$Nombre_cientifico <- NULL
  m <- t(as.matrix(abund))
  m <- m[, colnames(m) %in% species_ref, drop = FALSE]
  m <- m[, species_ref[species_ref %in% colnames(m)], drop = FALSE]
  empty <- rownames(m)[rowSums(m) == 0]
  if (length(empty) > 0) {
    log_msg(group, ": empty zones excluded -> ", paste(empty, collapse = ", "))
    m <- m[rowSums(m) > 0, , drop = FALSE]
  }
  m[, colSums(m) > 0, drop = FALSE]
}

# 4.3 FEASIBILITY DIAGNOSTIC ---------------------------------------------------
#  Which indices can each zone actually support? FRic, FEve and FDiv need the
#  convex hull, which needs S > number of PCoA axes retained (and at minimum
#  S >= 3). FDis does not. Publishing this table means an NA in the results is
#  a documented structural limit, not an unexplained gap -- which is exactly
#  what Zone 1 (2 plant species) needs.
feasibility_table <- function(comm, group, dir_out, n_axes_hint = NA) {
  S <- rowSums(comm > 0)
  tab <- data.frame(
    Zone      = rownames(comm),
    n_species = as.integer(S),
    FRic      = ifelse(S >= 3, "computable", "NA - convex hull needs S >= 3"),
    FEve      = ifelse(S >= 3, "computable", "NA - needs S >= 3"),
    FDiv      = ifelse(S >= 3, "computable", "NA - needs S >= 3"),
    FDis      = ifelse(S >= 1, "computable", "NA - empty zone"),
    note      = ifelse(S == 1, "single species: FDis = 0 by definition", ""),
    stringsAsFactors = FALSE)
  if (!is.na(n_axes_hint)) {
    tab$FRic <- ifelse(S > n_axes_hint, tab$FRic,
                       paste0("NA - S <= ", n_axes_hint, " retained PCoA axes"))
  }
  write.csv(tab, file.path(dir_out, paste0("index_feasibility_", group, ".csv")),
            row.names = FALSE)
  log_msg(group, " - richness per zone: ",
          paste(tab$Zone, tab$n_species, sep = "=", collapse = "  "))
  not_ok <- tab$Zone[tab$n_species < 3]
  if (length(not_ok) > 0)
    log_msg(group, " - zones where FRic/FEve/FDiv are NOT computable (S < 3): ",
            paste(not_ok, collapse = ", "))
  tab
}

# 4.4 Reusable dbFD core (used by the main analysis and by the sensitivity run).
dbFD_core <- function(traits, comm, group, tag = "") {
  res <- tryCatch(
    FD::dbFD(x = traits, a = comm, w.abun = TRUE, corr = pcoa_correction,
             calc.FRic = TRUE, calc.FDiv = TRUE, calc.CWM = FALSE,
             stand.x = TRUE, messages = FALSE),
    error = function(e) {
      log_msg(group, " [dbFD", tag, "] ERROR: ", conditionMessage(e)); NULL })
  if (is.null(res)) return(NULL)
  data.frame(Zone = rownames(comm), n_species = res$nbsp,
             FRic = res$FRic, FEve = res$FEve, FDiv = res$FDiv,
             FDis = res$FDis, RaoQ = res$RaoQ,
             qual_FRic = ifelse(is.null(res$qual.FRic), NA, res$qual.FRic))
}

# 4.5 ENGINE A — FD::dbFD -------------------------------------------------------
run_dbFD <- function(traits, comm, group, dir_out) {
  log_msg(group, " [dbFD] Gower + PCoA(", pcoa_correction, ")")
  d <- FD::gowdis(traits)
  euclid <- ade4::is.euclid(d)
  log_msg(group, " [dbFD] Gower Euclidean: ", euclid,
          if (!euclid) paste0(" -> ", pcoa_correction, " correction applied") else "")
  out <- dbFD_core(traits, comm, group)
  if (!is.null(out))
    write.csv(out, file.path(dir_out, paste0("indices_dbFD_", group, ".csv")),
              row.names = FALSE)
  out
}

# 4.6 ENGINE B — mFD (functional-space quality), optional ----------------------
run_mFD <- function(traits, comm, group, dir_out) {
  if (!mFD_available) {
    log_msg(group, " [mFD] skipped: package not installed (see section 1). ",
            "Analysis continues with FD::dbFD.")
    return(NULL)
  }
  log_msg(group, " [mFD] funct.dist + quality.fspaces + alpha.fd.multidim")
  types  <- sapply(traits, function(x) if (is.numeric(x)) "Q" else "N")
  tr_cat <- data.frame(trait_name = names(types), trait_type = as.character(types),
                       stringsAsFactors = FALSE)
  res <- tryCatch({
    sp_dist <- mFD::funct.dist(sp_tr = traits, tr_cat = tr_cat, metric = "gower")
    maxdim  <- min(nrow(traits) - 1, 6)
    qual    <- mFD::quality.fspaces(sp_dist = sp_dist, maxdim_pcoa = maxdim,
                                    deviation_weighting = "absolute",
                                    fdist_scaling = FALSE, fdendro = "average")
    qtab    <- qual$quality_fspaces
    mad_col <- qtab[[grep("mad|mAD", names(qtab), value = TRUE)[1]]]
    names(mad_col) <- rownames(qtab)
    mad_pcoa <- mad_col[grepl("^pcoa_", names(mad_col))]
    best     <- names(which.min(mad_pcoa))
    n_axes   <- as.integer(gsub("\\D", "", best))
    rich <- rowSums(comm > 0)
    if (n_axes >= max(rich)) n_axes <- max(1, max(rich) - 1)
    comm_mfd <- comm[rich > n_axes, , drop = FALSE]
    excl <- setdiff(rownames(comm), rownames(comm_mfd))
    if (length(excl) > 0)
      log_msg(group, " [mFD] zones excluded (richness <= ", n_axes,
              " axes): ", paste(excl, collapse = ", "))
    if (nrow(comm_mfd) == 0) stop("no zone with enough richness for mFD")
    comm_mfd <- comm_mfd[, colSums(comm_mfd) > 0, drop = FALSE]
    log_msg(group, " [mFD] chosen space: ", n_axes, " axes (min mAD=",
            round(min(mad_pcoa), 4), ")")
    coords <- qual$details_fspaces$sp_pc_coord[colnames(comm_mfd), 1:n_axes, drop = FALSE]
    alpha  <- mFD::alpha.fd.multidim(sp_faxes_coord = coords,
               asb_sp_w = as.matrix(comm_mfd),
               ind_vect = c("fric", "feve", "fdiv", "fdis"),
               scaling = TRUE, check_input = TRUE, details_returned = FALSE)
    list(idx = alpha$functional_diversity_indices, qtab = qtab)
  }, error = function(e) { log_msg(group, " [mFD] ERROR: ", conditionMessage(e)); NULL })
  if (is.null(res)) return(NULL)
  out <- res$idx; out$Zone <- rownames(out)
  write.csv(out, file.path(dir_out, paste0("indices_mFD_", group, ".csv")), row.names = FALSE)
  write.csv(res$qtab, file.path(dir_out, paste0("space_quality_mFD_", group, ".csv")))
  out
}

# 4.7 COMPARISON dbFD vs mFD ---------------------------------------------------
compare_engines <- function(db, mf, group, dir_out) {
  if (is.null(db) || is.null(mf)) return(invisible(NULL))
  m <- merge(
    data.frame(Zone = db$Zone, FRic_db = db$FRic, FEve_db = db$FEve,
               FDiv_db = db$FDiv, FDis_db = db$FDis),
    data.frame(Zone = mf$Zone, FRic_mFD = mf$fric, FEve_mFD = mf$feve,
               FDiv_mFD = mf$fdiv, FDis_mFD = mf$fdis),
    by = "Zone", all = TRUE)
  write.csv(m, file.path(dir_out, paste0("comparison_dbFD_vs_mFD_", group, ".csv")),
            row.names = FALSE)
  m
}

# 4.8 Label format: FRic in scientific notation (its values are tiny volumes
#     that would read as 0.000 at 3 decimals); the rest in plain decimal.
fmt_index <- function(value, index) {
  ifelse(is.na(value), "NA",
         ifelse(index == "FRic",
                formatC(value, format = "e", digits = 2),
                formatC(value, format = "f", digits = 3)))
}

# 4.9 HEAT MAP by Z-SCORE (relative contrasts between zones) -------------------
heatmap_zscore <- function(idx_df, cols_idx, file_name, dir_out, title,
                           row = "Zone") {
  df <- idx_df[, c(row, cols_idx)]
  z <- df
  for (cc in cols_idx) {
    v <- df[[cc]]
    s <- sd(v, na.rm = TRUE)                 # NA if fewer than 2 valid values
    z[[cc]] <- if (!is.na(s) && s > 0)
      (v - mean(v, na.rm = TRUE)) / s else rep(0, length(v))
  }
  lz <- reshape2::melt(z, id.vars = row, variable.name = "Index", value.name = "z")
  lv <- reshape2::melt(df, id.vars = row, variable.name = "Index", value.name = "value")
  la <- merge(lz, lv, by = c(row, "Index"))
  la$lab <- fmt_index(la$value, la$Index)
  g <- ggplot2::ggplot(la, ggplot2::aes(x = Index, y = .data[[row]], fill = z)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = lab), size = 3.1) +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                                  midpoint = 0, na.value = "grey85", name = "z-score") +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
  ggplot2::ggsave(file.path(dir_out, paste0(file_name, ".png")), g, width = 7, height = 4.4, dpi = 300)
  ggplot2::ggsave(file.path(dir_out, paste0(file_name, ".pdf")), g, width = 7, height = 4.4)
}

# 4.10 HEAT MAP by ABSOLUTE VALUE (real magnitudes) ----------------------------
heatmap_value <- function(idx_df, file_name, dir_out, title, row = "Zone") {
  panel <- function(cols, limits, pal, subtitle) {
    d <- reshape2::melt(idx_df[, c(row, cols), drop = FALSE], id.vars = row,
                        variable.name = "Index", value.name = "value")
    d$lab <- fmt_index(d$value, d$Index)
    ggplot2::ggplot(d, ggplot2::aes(x = Index, y = .data[[row]], fill = value)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.6) +
      ggplot2::geom_text(ggplot2::aes(label = lab), size = 3.1) +
      ggplot2::scale_fill_gradientn(colours = pal, limits = limits,
                                    na.value = "grey85", name = "value") +
      ggplot2::labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(panel.grid = ggplot2::element_blank())
  }
  max_fric <- suppressWarnings(max(idx_df$FRic, na.rm = TRUE))
  if (!is.finite(max_fric) || max_fric <= 0) max_fric <- 1
  g1 <- panel("FRic", c(0, max_fric), c("#F7FCF5", "#74C476", "#00441B"),
              "FRic (own scale: functional-space volume)")
  g2 <- panel(c("FEve", "FDiv", "FDis"), c(0, 1), c("#F7FBFF", "#6BAED6", "#08306B"),
              "FEve, FDiv, FDis (0-1 scale)")
  ggplot2::ggsave(file.path(dir_out, paste0(file_name, "_FRic.png")), g1, width = 3.6, height = 4.4, dpi = 300)
  ggplot2::ggsave(file.path(dir_out, paste0(file_name, "_FEve_FDiv_FDis.png")), g2, width = 6, height = 4.4, dpi = 300)
  ggplot2::ggsave(file.path(dir_out, paste0(file_name, "_FRic.pdf")), g1, width = 3.6, height = 4.4)
  ggplot2::ggsave(file.path(dir_out, paste0(file_name, "_FEve_FDiv_FDis.pdf")), g2, width = 6, height = 4.4)
}

# 4.11 SES heat map ------------------------------------------------------------
heatmap_ses <- function(ses_df, file_name, dir_out, title) {
  d <- ses_df
  d$lab <- ifelse(is.na(d$SES), "NA",
                  paste0(formatC(d$SES, format = "f", digits = 2),
                         ifelse(!is.na(d$p_two_tailed) & d$p_two_tailed < 0.05, "*", "")))
  g <- ggplot2::ggplot(d, ggplot2::aes(x = Index, y = Zone, fill = SES)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = lab), size = 3.1) +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                                  midpoint = 0, na.value = "grey85", name = "SES") +
    ggplot2::labs(title = title,
                  subtitle = "SES vs 999 trait-label permutations; * = p < 0.05.\nPositive = overdispersed, negative = clustered.",
                  x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
  ggplot2::ggsave(file.path(dir_out, paste0(file_name, ".png")), g, width = 7.5, height = 4.6, dpi = 300)
  ggplot2::ggsave(file.path(dir_out, paste0(file_name, ".pdf")), g, width = 7.5, height = 4.6)
}

# 4.12 ABUNDANCE-WEIGHTING SENSITIVITY -----------------------------------------
#  Recompute indices with raw / log1p / presence-absence weights. FRic is
#  identical across the three (it depends only on presence); large gaps between
#  'raw' and 'presence' for FEve/FDiv/FDis flag dominance by one or few species.
sensitivity_abundance <- function(traits, comm, group, dir_out) {
  dir_sens <- file.path(dir_out, "sensitivity")
  dir.create(dir_sens, recursive = TRUE, showWarnings = FALSE)
  log_msg(group, " [sensitivity] raw vs log1p vs presence")
  weightings <- list(raw = comm, log1p = log1p(comm), presence = (comm > 0) * 1)
  tabs <- list()
  for (nm in names(weightings)) {
    out <- dbFD_core(traits, weightings[[nm]], group, paste0("/", nm))
    if (!is.null(out)) { out$Weighting <- nm; tabs[[nm]] <- out }
  }
  if (length(tabs) == 0) return(invisible(NULL))
  tab <- do.call(rbind, tabs)
  write.csv(tab, file.path(dir_sens, paste0("sensitivity_", group, ".csv")), row.names = FALSE)
  base_r <- tab[tab$Weighting == "raw",      c("Zone", "FDis")]
  base_p <- tab[tab$Weighting == "presence", c("Zone", "FDis")]
  delta  <- merge(base_r, base_p, by = "Zone", suffixes = c("_raw", "_presence"))
  delta$delta_FDis <- delta$FDis_presence - delta$FDis_raw
  write.csv(delta, file.path(dir_sens, paste0("delta_FDis_", group, ".csv")), row.names = FALSE)
  dominated <- delta$Zone[abs(delta$delta_FDis) > 0.15 & !is.na(delta$delta_FDis)]
  if (length(dominated) > 0)
    log_msg(group, " [sensitivity] zones with strong dominance effect ",
            "(|delta FDis| > 0.15): ", paste(dominated, collapse = ", "))
  tab$ZW <- paste(tab$Zone, tab$Weighting, sep = " . ")
  heatmap_zscore(tab[, c("ZW","FRic","FEve","FDiv","FDis")],
                 c("FRic","FEve","FDiv","FDis"),
                 paste0("heatmap_sensitivity_", group), dir_sens,
                 paste0("Weighting sensitivity - ", group), row = "ZW")
  invisible(tab)
}

# 4.13 ORCHESTRATOR: run everything for one group ------------------------------
analyze_group <- function(traits, comm, group, dir_out, title, counts = TRUE) {
  log_msg("========== ", group, " ==========")
  feasibility_table(comm, group, dir_out)
  write.csv(data.frame(Species = rownames(traits), traits),
            file.path(dir_out, paste0("trait_matrix_", group, ".csv")), row.names = FALSE)
  write.csv(data.frame(Zone = rownames(comm), comm),
            file.path(dir_out, paste0("community_matrix_", group, ".csv")), row.names = FALSE)

  db <- run_dbFD(traits, comm, group, dir_out)
  mf <- run_mFD (traits, comm, group, dir_out)
  compare_engines(db, mf, group, dir_out)

  # --- functional redundancy (Ricotta et al. 2016) ---------------------------
  d_gow <- FD::gowdis(traits)
  red   <- functional_redundancy(comm, d_gow)
  write.csv(red, file.path(dir_out, paste0("redundancy_", group, ".csv")), row.names = FALSE)
  log_msg(group, " [redundancy] relative FR per zone: ",
          paste(red$Zone, formatC(red$FR_rel, format = "f", digits = 3),
                sep = "=", collapse = "  "))

  # --- taxonomic diversity, and the taxonomic-vs-functional contrast ---------
  tax <- taxonomic_diversity(comm, counts = counts)
  write.csv(tax, file.path(dir_out, paste0("taxonomic_diversity_", group, ".csv")),
            row.names = FALSE)
  if (!is.null(db)) {
    tvf <- taxonomic_vs_functional(tax, db[, c("Zone","FRic","FEve","FDiv","FDis")], red)
    write.csv(tvf, file.path(dir_out, paste0("taxonomic_vs_functional_", group, ".csv")),
              row.names = FALSE)
  }

  # --- null models -> richness-independent effect sizes ----------------------
  ses <- NULL
  if (run_null_models) {
    log_msg(group, " [null models] ", n_perm, " trait-label permutations ",
            "(this is the slow step)")
    ses <- tryCatch(
      ses_functional(traits, comm, n_perm = n_perm, corr = pcoa_correction),
      error = function(e) { log_msg(group, " [null models] ERROR: ",
                                    conditionMessage(e)); NULL })
    if (!is.null(ses)) {
      write.csv(ses, file.path(dir_out, paste0("ses_null_models_", group, ".csv")),
                row.names = FALSE)
      heatmap_ses(ses, paste0("heatmap_SES_", group), dir_out,
                  paste0("SES vs null model - ", group))
      sig <- ses[!is.na(ses$p_two_tailed) & ses$p_two_tailed < 0.05, ]
      if (nrow(sig) > 0) {
        log_msg(group, " [null models] departures from random: ",
                paste(sig$Zone, sig$Index, sep = "/", collapse = ", "))
      } else {
        log_msg(group, " [null models] no index departs from the null in any ",
                "zone: differences between zones are consistent with richness ",
                "differences alone.")
      }
    }
  }

  # --- heat maps of the raw indices ------------------------------------------
  if (!is.null(db)) {
    heatmap_zscore(db, c("FRic","FEve","FDiv","FDis"), paste0("heatmap_", group),
                   dir_out, title)
    heatmap_value(db, paste0("heatmap_value_", group), dir_out, title)
  } else if (!is.null(mf)) {
    mf2 <- data.frame(Zone = mf$Zone, FRic = mf$fric, FEve = mf$feve,
                      FDiv = mf$fdiv, FDis = mf$fdis)
    heatmap_zscore(mf2, c("FRic","FEve","FDiv","FDis"), paste0("heatmap_", group),
                   dir_out, title)
    heatmap_value(mf2, paste0("heatmap_value_", group), dir_out, title)
  }

  sensitivity_abundance(traits, comm, group, dir_out)
  invisible(list(dbFD = db, mFD = mf, redundancy = red, taxonomic = tax, ses = ses))
}

## ---- 5. RUN ------------------------------------------------------------------

# 5.1 FLORA - 7 traits: H (log10), LDMC, SLA (log10) + 4 categorical -----------
#     counts = FALSE: flora abundance is crown cover, so rarefaction is not
#     meaningful and is deliberately not computed.
R_flora <- prep_traits(traits_flora,
             continuous  = c("H_m", "LDMC", "SLA"),
             categorical = c("forma_crecimiento","fenologia","forma_hoja","dispersion"),
             log_cols    = if (use_log_height) c("H_m", "SLA") else "SLA")
C_flora <- prep_comm(abund_flora, rownames(R_flora), "flora")
R_flora <- R_flora[colnames(C_flora), , drop = FALSE]
res_flora <- analyze_group(R_flora, C_flora, "flora", dir_flora,
               "Functional diversity of FLORA - Las Penas Campus", counts = FALSE)

# 5.2 FAUNA - all groups, three campaigns --------------------------------------
#     Kept because it is the analysis that MOTIVATES the birds-only run below:
#     the Apis mellifera dominance shows up here as a depressed FDis in Zone 4.
fauna_list <- list(
  list(traits = traits_fau[[1]], abund = abund_fau[[1]], id = "fauna_c1",
       dir = dirs_fau[1], title = "Functional diversity of FAUNA - Campaign 1"),
  list(traits = traits_fau[[2]], abund = abund_fau[[2]], id = "fauna_c2",
       dir = dirs_fau[2], title = "Functional diversity of FAUNA - Campaign 2"),
  list(traits = traits_fau[[3]], abund = abund_fau[[3]], id = "fauna_c3",
       dir = dirs_fau[3], title = "Functional diversity of FAUNA - Campaign 3"))

res_fauna <- list()
for (cf in fauna_list) {
  R_f <- prep_traits(cf$traits, continuous = "masa_g",
                     categorical = c("gremio_trofico","estrato"), log_cols = "masa_g")
  C_f <- prep_comm(cf$abund, rownames(R_f), cf$id)
  R_f <- R_f[colnames(C_f), , drop = FALSE]
  res_fauna[[cf$id]] <- analyze_group(R_f, C_f, cf$id, cf$dir, cf$title)
}

# 5.3 BIRDS ONLY, CAMPAIGNS POOLED ---------------------------------------------
#  This is the analysis the thesis reports (Table 9), and the reason for it is
#  substantive, not cosmetic. Apis mellifera accounted for 90.1-95.7 % of the
#  Zone 4 fauna records, and that figure is an EXTRAPOLATED estimate on a
#  different scale from the vertebrate detections (see data/README.md). Because
#  FDis weights distances by relative abundance, pooling it with the birds drags
#  the community centroid onto a single functional type and depresses the index
#  -- a methodological artefact, not a property of the community.
#  Campaigns are pooled because per campaign and zone the richness is 1-4
#  species, below what any of these indices can support.
birds <- sp_meta$Nombre_cientifico[sp_meta$grupo == "Ave"]
log_msg("========== birds_pooled: ", length(birds), " bird species ==========")

abund_birds <- Reduce(function(a, b) {
  m <- merge(a, b, by = "Nombre_cientifico", all = TRUE)
  zones <- unique(sub("\\.[xy]$", "", setdiff(names(m), "Nombre_cientifico")))
  out <- data.frame(Nombre_cientifico = m$Nombre_cientifico)
  for (z in zones) {
    cols <- grep(paste0("^", z, "(\\.[xy])?$"), names(m), value = TRUE)
    out[[z]] <- rowSums(m[, cols, drop = FALSE], na.rm = TRUE)
  }
  out
}, lapply(abund_fau, function(a) a[a$Nombre_cientifico %in% birds, , drop = FALSE]))

traits_birds <- do.call(rbind, traits_fau)
traits_birds <- traits_birds[!duplicated(traits_birds$Nombre_cientifico), ]
traits_birds <- traits_birds[traits_birds$Nombre_cientifico %in% birds, ]

R_b <- prep_traits(traits_birds, continuous = "masa_g",
                   categorical = c("gremio_trofico","estrato"), log_cols = "masa_g")
C_b <- prep_comm(abund_birds, rownames(R_b), "birds_pooled")
R_b <- R_b[colnames(C_b), , drop = FALSE]
res_birds <- analyze_group(R_b, C_b, "birds_pooled", dir_birds,
               "Functional diversity of BIRDS - campaigns pooled")

# 5.4 CROSS-CAMPAIGN FAUNA COMPARISON ------------------------------------------
tabs <- lapply(names(res_fauna), function(id) {
  db <- res_fauna[[id]]$dbFD
  if (is.null(db)) return(NULL)
  db$Campaign <- gsub("fauna_", "C", id)   # C1, C2, C3
  db
})
comp <- do.call(rbind, tabs[!sapply(tabs, is.null)])
if (!is.null(comp)) {
  write.csv(comp, file.path(dir_fcomp, "fauna_indices_by_campaign.csv"), row.names = FALSE)
  comp$SZ <- paste(comp$Campaign, comp$Zone, sep = " . ")
  heatmap_zscore(comp[, c("SZ","FRic","FEve","FDiv","FDis")],
                 c("FRic","FEve","FDiv","FDis"), "heatmap_fauna_campaigns",
                 dir_fcomp, "FAUNA - indices by campaign and zone", row = "SZ")
  log_msg("Cross-campaign comparison saved to ", dir_fcomp)
}

# 5.5 APIS MELLIFERA DOMINANCE DIAGNOSTIC --------------------------------------
#  Quantify, per campaign and zone, the share of fauna records held by the most
#  abundant species. This is the number behind the decision to analyse birds
#  separately, so it belongs in the output rather than only in the discussion.
apis_rows <- list()
for (i in 1:3) {
  a <- abund_fau[[i]]
  zones <- setdiff(names(a), "Nombre_cientifico")
  for (z in zones) {
    v <- a[[z]]; tot <- sum(v)
    if (tot == 0) next
    k <- which.max(v)
    apis_rows[[length(apis_rows) + 1L]] <- data.frame(
      Campaign = paste0("C", i), Zone = z, total_records = tot,
      dominant_species = a$Nombre_cientifico[k],
      dominant_share_pct = round(100 * max(v) / tot, 1),
      apis_share_pct = round(100 * sum(v[a$Nombre_cientifico == "Apis mellifera"]) / tot, 1))
  }
}
dom <- do.call(rbind, apis_rows)
write.csv(dom, file.path(dir_base, "fauna", "dominance_diagnostic.csv"), row.names = FALSE)
log_msg("Dominance diagnostic: max Apis share = ",
        max(dom$apis_share_pct, na.rm = TRUE), " % (Zone 4)")

## ---- 6. WRAP-UP --------------------------------------------------------------
writeLines(capture.output(sessionInfo()), file.path(dir_logs, "sessionInfo.txt"))
log_msg("Analysis finished. Results in: ", dir_base)
cat("\n>>> DONE. Check the folder:", dir_base, "\n")
