# ==============================================================================
#  NULL MODELS — richness-independent standardized effect sizes (SES)
# ------------------------------------------------------------------------------
#  THE PROBLEM THIS SOLVES.
#  FRic, FDis and Rao's Q all increase with species richness by construction.
#  In this study the assessed zones range from 2 to 26 plant species, so a raw
#  comparison between zones confounds "how many species are there" with "how
#  the species are arranged in trait space". Saying "Zone 4 is functionally
#  richer than Zone 3" is, without a null model, close to saying "Zone 4 has
#  more species" -- which we already knew from the census.
#
#  THE FIX.
#  Compare each observed index against a null distribution built by shuffling
#  the species labels of the TRAIT matrix while holding the community matrix
#  (and therefore richness and the abundance distribution) fixed. This is the
#  classic "taxa.labels" / independent-swap-on-traits null of Gotelli & Graves
#  (1996) and Swenson (2014): it asks whether the species that co-occur in a
#  zone are more (or less) functionally similar than an equally rich, equally
#  abundance-structured random draw from the campus species pool.
#
#       SES = (observed - mean(null)) / sd(null)
#
#  SES > 0  -> functional OVERDISPERSION (co-occurring species are more
#              different than expected). Often read as limiting similarity /
#              niche differentiation.
#  SES < 0  -> functional CLUSTERING (co-occurring species are more similar
#              than expected). Often read as environmental filtering -- the
#              expected signature of urbanization.
#  |SES| > ~1.96 is the conventional two-tailed 5 % threshold, but the
#  permutational p-value reported here is exact and should be preferred.
#
#  WHAT THE POOL IS.
#  The pool is the set of species with traits in the group being analysed
#  (all 42 plants, or all birds). It is a CAMPUS-level pool, not a regional
#  one, so the result answers "given the species present on this campus, is
#  this zone's assemblage a random subset?" -- not "is the campus itself
#  filtered relative to the dry forest?". The second question needs a regional
#  species pool that this study did not sample, and the limitation is stated in
#  the README rather than hidden in the code.
# ==============================================================================

#' Standardized effect sizes for functional diversity indices
#'
#' @param traits trait data.frame, rownames = species (mixed types allowed)
#' @param comm community matrix, sites x species
#' @param n_perm number of permutations (999 gives p-values to 0.001)
#' @param corr PCoA correction passed to FD::dbFD
#' @param indices character vector of index names to standardize
#' @param seed integer seed, set inside the function so a single call is
#'   reproducible on its own
#' @param quiet suppress dbFD's own messages
#' @return data.frame, one row per site x index, with observed value, null
#'   mean and sd, SES, and a two-tailed permutational p-value
ses_functional <- function(traits, comm, n_perm = 999,
                           corr = "cailliez",
                           indices = c("FRic", "FEve", "FDiv", "FDis", "RaoQ"),
                           seed = 123, quiet = TRUE) {

  if (!requireNamespace("FD", quietly = TRUE))
    stop("Package 'FD' is required for ses_functional().")

  set.seed(seed)

  run_once <- function(tr) {
    res <- try(FD::dbFD(x = tr, a = comm, w.abun = TRUE, corr = corr,
                        calc.FRic = TRUE, calc.FDiv = TRUE, calc.CWM = FALSE,
                        stand.x = TRUE, messages = FALSE), silent = TRUE)
    if (inherits(res, "try-error")) return(NULL)
    vapply(indices, function(i) {
      v <- res[[i]]
      if (is.null(v)) rep(NA_real_, nrow(comm)) else as.numeric(v)
    }, numeric(nrow(comm)))
  }

  obs <- run_once(traits)
  if (is.null(obs))
    stop("ses_functional(): the observed run of dbFD failed; ",
         "nothing to standardize against.")

  # Null distribution: shuffle species labels on the TRAIT matrix only.
  # Richness, abundances and the community matrix are untouched by construction.
  null_arr <- array(NA_real_, dim = c(nrow(comm), length(indices), n_perm))
  n_ok <- 0L
  for (k in seq_len(n_perm)) {
    tr_perm <- traits
    rownames(tr_perm) <- sample(rownames(traits))
    tr_perm <- tr_perm[colnames(comm), , drop = FALSE]
    r <- run_once(tr_perm)
    if (!is.null(r)) { null_arr[, , k] <- r; n_ok <- n_ok + 1L }
  }
  if (n_ok < n_perm * 0.5)
    warning("ses_functional(): only ", n_ok, "/", n_perm,
            " permutations succeeded; SES values are unreliable.")

  rows <- list()
  for (s in seq_len(nrow(comm))) {
    for (j in seq_along(indices)) {
      o  <- obs[s, j]
      nl <- null_arr[s, j, ]
      nl <- nl[is.finite(nl)]
      if (!is.finite(o) || length(nl) < 2) {
        rows[[length(rows) + 1L]] <- data.frame(
          Zone = rownames(comm)[s], Index = indices[j], observed = o,
          null_mean = NA_real_, null_sd = NA_real_, SES = NA_real_,
          p_two_tailed = NA_real_, n_perm_valid = length(nl))
        next
      }
      m <- mean(nl); s_ <- stats::sd(nl)
      # Two-tailed exact permutational p, with the observed value counted in
      # both tails (Davison & Hinkley 1997) so p is never 0.
      p_lo <- (sum(nl <= o) + 1) / (length(nl) + 1)
      p_hi <- (sum(nl >= o) + 1) / (length(nl) + 1)
      rows[[length(rows) + 1L]] <- data.frame(
        Zone = rownames(comm)[s], Index = indices[j], observed = o,
        null_mean = m, null_sd = s_,
        SES = if (is.finite(s_) && s_ > 0) (o - m) / s_ else NA_real_,
        p_two_tailed = min(1, 2 * min(p_lo, p_hi)),
        n_perm_valid = length(nl))
    }
  }
  out <- do.call(rbind, rows)
  out$interpretation <- ifelse(
    is.na(out$SES), "not computable",
    ifelse(out$p_two_tailed >= 0.05, "indistinguishable from random",
           ifelse(out$SES > 0, "overdispersed (more different than expected)",
                  "clustered (more similar than expected)")))
  rownames(out) <- NULL
  out
}
