# ─────────────────────────────────────────────────────────────
#  FUNCTION TO CALCULATE THE “ARCTIC N-P” TRACER (Newton et al., 2013)
# ─────────────────────────────────────────────────────────────
# This function estimates the relative influence of Atlantic vs. Pacific Water
# using the orthogonal distance of a data point (N, P) to the Atlantic and Pacific
# end-member N:P lines. Returns a normalized value:
#   ANP = 0 → pure Atlantic end-member
#   ANP = 1 → pure Pacific end-member
ANP <- function(Phosphate, Nitrogen,
                slope_AW  = 18.34862,  int_AW  = -3.513761, # Atlantic line (default: Jones et al.)
                slope_PW  = 13.957,  int_PW  = -11.306) # Pacific line (default: Yamamoto-Kawai et al.)
{
  stopifnot(length(Phosphate) == length(Nitrogen)) # Ensure vectors are the same length

  # Compute the orthogonal distance from a point (x, y) to a line y = a*x + b
  d_line <- function(x, y, a, b) abs(a * x - y + b) / sqrt(a^2 + 1)

  d_AW <- d_line(Phosphate, Nitrogen, slope_AW, int_AW) # Distance to Atlantic line
  d_PW <- d_line(Phosphate, Nitrogen, slope_PW, int_PW) # Distance to Pacific line

  ANP <- d_AW / (d_AW + d_PW) # Normalized ratio
  return(ANP)
}

# ─────────────────────────────────────────────────────────────
#  FUNCTION TO CALCULATE PACIFIC WATER FRACTION (fPW)
# ─────────────────────────────────────────────────────────────
# This function estimates the fraction of Pacific Water (fPW) based on measured
# phosphate and nitrogen, using user-specified or default end-member lines.
#   fPW = 0 → pure Atlantic
#   fPW = 1 → pure Pacific
# The formula uses the relative position between the Atlantic and Pacific lines.
fpw <- function(Phosphate, Nitrogen,
                df_ref = NULL, # Optional: table of line coefficients
                atlantic_ref, # e.g. "Jones1998"
                pacific_ref, # e.g. "Yamamoto.Kawai2008"
                clip = TRUE) {
  # Provide default reference lines if none supplied
  if (is.null(df_ref)) {
    df_ref <- tibble::tribble(
      ~Water, ~Ref, ~slope, ~intercept,
      # "Atlantic", "Jones1998", 17.499, -3.07,
      "Atlantic", "Sutherland2009", 20.7, -2.70,
      "Atlantic", "Sherwood2021", 15.54, -0.26,
      # "Atlantic", "Dodd2012", 18.87, -3.21,
      "Atlantic", "Jones2008", 18.34862, -3.513761,
      # "Pacific", "Jones1998", 12.368, -10.54,
      "Pacific", "Sutherland2009", 12.4, -10.5,
      "Pacific", "Yamamoto.Kawai2008", 13.957, -11.306,
      "Pacific", "Sherwood2021", 14.07, -11.53,
      # "Pacific", "Dodd2012", 15.38, -14.46,
      "Pacific", "Jones2008", 15.31394, -14.3951
    )
  }

  # Get line coefficients for each end-member
  atl_row <- df_ref %>% filter(Water == "Atlantic", Ref == atlantic_ref)
  pac_row <- df_ref %>% filter(Water == "Pacific", Ref == pacific_ref)
  stopifnot(nrow(atl_row) == 1, nrow(pac_row) == 1)

  slope_AW <- atl_row$slope
  int_AW <- atl_row$intercept
  slope_PW <- pac_row$slope
  int_PW <- pac_row$intercept

  # Compute N for each end-member at measured phosphate
  Nitrogen_AW_star <- slope_AW * Phosphate + int_AW
  Nitrogen_PW_star <- slope_PW * Phosphate + int_PW

  # Compute fPW as the fractional distance between Atlantic and Pacific lines
  fpw_raw <- (Nitrogen - Nitrogen_AW_star) / (Nitrogen_PW_star - Nitrogen_AW_star)
  if (clip) fpw_raw <- pmin(pmax(fpw_raw, 0), 1) # Optional: force result between 0 and 1

  return(fpw_raw)
}

# ─────────────────────────────────────────────────────────────
#  FUNCTION TO CALCULATE THE POs* TRACER (Jones et al., 2008)
# ─────────────────────────────────────────────────────────────
# This function calculates the phosphate–oxygen (POs*) tracer, which integrates
# phosphate, oxygen, and salinity to distinguish water masses in the Arctic.
# The tracer transitions from a Redfield O:P ratio above the halocline (175:1)
# to a lower value below (125:1), based on salinity.
#
# Arguments:
#   phosphate - vector of phosphate concentration (µmol/kg)
#   oxygen    - vector of oxygen concentration (µmol/kg)
#   salinity  - vector of practical salinity (unitless)
# Returns:
#   Vector of POs* values
POs_star <- function(phosphate, oxygen, salinity) {
  # S-shaped weighting: transition at S = 34.5
  f <- 1 / (1 + exp(10 * (34.5 - salinity)))

  # Above halocline (Redfield: 175:1)
  POs1 <- phosphate + oxygen / 175 - 1.95
  # Below halocline (125:1)
  POs2 <- phosphate + oxygen / 125 - 1.95

  # Weighted mean for a smooth transition
  POs_star <- f * POs1 + (1 - f) * POs2
  return(POs_star)
}

# Example usage for POs*:
# df$POs_star <- POs_star(df$Phosphate, df$Oxygen, df$Salinity)


