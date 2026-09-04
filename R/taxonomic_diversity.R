# ==============================================================================
#  TAXONOMIC DIVERSITY — the contrast that makes the functional result readable
# ------------------------------------------------------------------------------
#  WHY THIS IS HERE.
#  "Functional richness is low" only becomes an argument once you can say
#  "...even though species richness is not". The whole resilience claim of this
#  study rests on that gap: many species, few distinct functions, therefore high
#  redundancy and some buffering capacity against local species loss.
#  Reporting taxonomic and functional diversity side by side is what turns four
#  index values into that statement.
#
#  WHAT IS COMPUTED.
#    S           species richness
#    Shannon H'  -sum(p * log(p))            (sensitive to rare species)
#    Simpson D   1 - sum(p^2)                (sensitive to dominant species)
#    Pielou J'   H' / log(S)                 (evenness, 0-1)
#    Berger-Parker  max(p)                   (plain dominance; the Apis diagnostic)
#    Rarefied S  expected richness at a common sample size (COUNTS ONLY)
#
#  ON RAREFACTION. It is only meaningful for integer counts of individuals.
#  Flora abundance here is CROWN COVER, a continuous proportion, so rarefaction
#  is not computed for flora and the column is returned as NA rather than
#  silently coercing cover to counts. For fauna, note that the sampling effort
#  per unit area varies 4.4-fold between zones (see data/zones.csv), so
#  rarefaction to the smallest zone total is the honest comparison, and even
#  that does not fully correct for the effort difference.
# ==============================================================================

#' Taxonomic diversity indices per site
#'
#' @param comm community matrix, sites x species
#' @param counts logical; TRUE if the values are integer counts of individuals
#'   (enables rarefaction), FALSE for cover or other continuous abundance
#' @return data.frame with one row per site
taxonomic_diversity <- function(comm, counts = TRUE) {

  rare_to <- NA_real_
  if (counts) {
    tot <- rowSums(comm)
    tot <- tot[tot > 0]
    if (length(tot) > 0) rare_to <- min(tot)
  }

  rows <- lapply(rownames(comm), function(site) {
    a <- as.numeric(comm[site, ])
    a <- a[a > 0]
    S <- length(a)
    if (S == 0) {
      return(data.frame(Zone = site, S = 0L, Shannon = NA_real_,
                        Simpson = NA_real_, Pielou = NA_real_,
                        BergerParker = NA_real_, S_rarefied = NA_real_))
    }
    p <- a / sum(a)
    H <- -sum(p * log(p))
    Sr <- NA_real_
    if (counts && is.finite(rare_to) &&
        requireNamespace("vegan", quietly = TRUE) &&
        all(abs(a - round(a)) < 1e-8)) {
      Sr <- as.numeric(vegan::rarefy(round(a), sample = rare_to))
    }
    data.frame(Zone         = site,
               S            = S,
               Shannon      = H,
               Simpson      = 1 - sum(p^2),
               # Pielou is undefined for a single species (log(1) = 0):
               # report NA, not 1, because there is no evenness to speak of.
               Pielou       = if (S > 1) H / log(S) else NA_real_,
               BergerParker = max(p),
               S_rarefied   = Sr)
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  attr(out, "rarefied_to") <- rare_to
  out
}

#' Join taxonomic and functional diversity into one interpretive table
#'
#' The point of this table is the last column: how much of the species
#' diversity is NOT expressed as functional diversity.
#'
#' @param tax output of taxonomic_diversity()
#' @param fd  data.frame with a Zone column and FRic/FEve/FDiv/FDis
#' @param red output of functional_redundancy() (optional)
taxonomic_vs_functional <- function(tax, fd, red = NULL) {
  out <- merge(tax, fd, by = "Zone", all = TRUE)
  if (!is.null(red)) {
    out <- merge(out, red[, c("Zone", "RaoQ", "FR_rel")], by = "Zone", all = TRUE)
    out$pct_diversity_not_functional <- round(100 * out$FR_rel, 1)
  }
  out[order(out$Zone), ]
}
