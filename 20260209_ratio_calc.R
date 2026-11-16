# ------------------------------------------------------------------------------
# safe_ratio.R
#
# Robust computation of nutrient stoichiometric ratios (N:P, Si:N, Si:P)
# for the NOW transect (Baffin Bay / Nares Strait).
#
# Nutrient ratios are undefined when either the numerator or the denominator
# falls below its analytical limit of detection (LOD). In such cases the ratio
# is returned as NA to avoid spurious infinite (Inf) or zero values that would
# otherwise bias water-mass-level summary statistics.
#
# Author : Nicolas Schiffrine
# Date   : 2026-06-18
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Analytical limits of detection (LOD), in micromoles per litre (umol L-1).
# Update these values to match your laboratory's certified detection limits.
# ------------------------------------------------------------------------------
LOD <- c(
  Nitrate   = 0.05,
  Phosphate = 0.05,
  Silicate  = 0.10
)

# ------------------------------------------------------------------------------
#' Compute a detection-limit-aware nutrient ratio
#'
#' Returns num / den, but sets the result to NA whenever the numerator or the
#' denominator falls below its respective limit of detection. This prevents
#' undefined (Inf) or non-informative ratios arising from near-zero nutrient
#' concentrations (e.g. summer surface nitrate drawdown).
#'
#' @param num     Numeric vector. Numerator nutrient concentration (umol L-1).
#' @param den     Numeric vector. Denominator nutrient concentration (umol L-1).
#' @param lod_num Numeric scalar. Limit of detection for the numerator nutrient.
#' @param lod_den Numeric scalar. Limit of detection for the denominator nutrient.
#'
#' @return Numeric vector of ratios, with NA where num < lod_num or den < lod_den.
#'
#' @examples
#' safe_ratio(c(5, 0, 10), c(0.5, 0.8, 0), LOD["Nitrate"], LOD["Phosphate"])
# ------------------------------------------------------------------------------
ratio_calc <- function(num, den, lod_num, lod_den) {
  dplyr::if_else(num < lod_num | den < lod_den, NA_real_, num / den)
}

