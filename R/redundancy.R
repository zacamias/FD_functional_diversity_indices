# ==============================================================================
#  FUNCTIONAL REDUNDANCY  (Ricotta et al. 2016, Methods Ecol. Evol. 7:1386-1395)
# ------------------------------------------------------------------------------
#  Functional redundancy is the fraction of SPECIES diversity that does NOT
#  translate into FUNCTIONAL diversity:
#
#       FR      = D - Q          (absolute redundancy)
#       FR_rel  = 1 - Q / D      (relative redundancy, 0-1)
#
#  where  D = Simpson diversity          = 1 - sum(p_i^2)
#         Q = Rao's quadratic entropy    = sum_i sum_j d_ij * p_i * p_j
#
#  and d_ij is a dissimilarity in [0, 1]. Gower distance already is; any other
#  distance is rescaled by its maximum before use, and the rescaling is
#  reported so the value is never silently altered.
#
#  WHY THIS MATTERS HERE. A low FRic tells you the community occupies a small
#  volume of trait space. It does NOT tell you that different species do the
#  same job -- that is a separate claim, and it needs its own measurement.
#  Redundancy is the quantity that supports the resilience argument: a
#  community where many species are functionally interchangeable can absorb the
#  loss of some of them without losing the function.
#
#  Q = 0 (all species functionally identical) gives FR_rel = 1.
#  Q = D (every pair maximally distinct) gives FR_rel = 0.
# ==============================================================================

#' Rao's quadratic entropy for one community
#'
#' @param p numeric vector of relative abundances (must sum to 1)
#' @param d square dissimilarity matrix, values in [0, 1], aligned to `p`
rao_Q <- function(p, d) {
  as.numeric(t(p) %*% d %*% p)
}

#' Simpson diversity (Gini-Simpson index) for one community
simpson_D <- function(p) {
  1 - sum(p^2)
}

#' Functional redundancy per community
#'
#' @param comm community matrix, sites (rows) x species (columns)
#' @param dist_sp species dissimilarity ("dist" object or matrix), species names
#'   must match `colnames(comm)`
#' @return data.frame with one row per site:
#'   \describe{
#'     \item{Simpson_D}{Species diversity (Gini-Simpson).}
#'     \item{RaoQ}{Rao's quadratic entropy computed on the rescaled distance.}
#'     \item{FR_abs}{D - Q. Interpretable only against the D of the same site.}
#'     \item{FR_rel}{1 - Q/D. Comparable across sites. NA when D = 0
#'       (monospecific community: redundancy is undefined, not zero).}
#'     \item{n_species}{Observed species richness.}
#'   }
functional_redundancy <- function(comm, dist_sp) {
  d <- as.matrix(dist_sp)

  # Align distance matrix to the community matrix, dropping species absent from
  # the trait set rather than silently recycling.
  sp <- intersect(colnames(comm), colnames(d))
  if (length(sp) < ncol(comm)) {
    warning("functional_redundancy: ", ncol(comm) - length(sp),
            " species in the community matrix have no trait distance; dropped.")
  }
  d    <- d[sp, sp, drop = FALSE]
  comm <- comm[, sp, drop = FALSE]

  # Rescale to [0, 1] if needed. Gower is already bounded; a corrected PCoA
  # distance may not be.
  dmax <- max(d, na.rm = TRUE)
  if (is.finite(dmax) && dmax > 1) {
    d <- d / dmax
    attr(d, "rescaled_by") <- dmax
  }

  out <- lapply(rownames(comm), function(site) {
    a <- comm[site, ]
    tot <- sum(a)
    if (tot <= 0) {
      return(data.frame(Zone = site, n_species = 0L, Simpson_D = NA_real_,
                        RaoQ = NA_real_, FR_abs = NA_real_, FR_rel = NA_real_))
    }
    p <- a / tot
    D <- simpson_D(p)
    Q <- rao_Q(p, d)
    data.frame(Zone      = site,
               n_species = sum(a > 0),
               Simpson_D = D,
               RaoQ      = Q,
               FR_abs    = D - Q,
               # D = 0 means a single species holds all the abundance; the
               # ratio is undefined and must not be reported as 0 or 1.
               FR_rel    = if (D > 0) 1 - Q / D else NA_real_)
  })

  res <- do.call(rbind, out)
  rownames(res) <- NULL
  attr(res, "distance_rescaled_by") <- if (is.finite(dmax) && dmax > 1) dmax else 1
  res
}
