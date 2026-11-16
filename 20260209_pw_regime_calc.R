#' Calculate Zonal Pacific Water Indices with a 6-Tier Regime Classification
#'
#' Analyses the spatiotemporal distribution of Pacific Water (PW) in the North
#' Water (NOW) polynya from the Pacific Water fraction (fPW) and classifies each
#' year into one of six hydrographic regimes.
#'
#' @param data Dataframe with columns: classification, Depth, fPW (0-1),
#'   Station_Cluster, decimalLongitude (negative = West), year.
#' @param min_stations Minimum number of unique stations required (e.g., 6).
#' @param fPW_cutoff fPW threshold defining a PW-dominated sample (default 0.50).
#' @param front_detection_thr fPW threshold defining the PW leading edge
#'   (default 0.40).
#' @param western_limit Longitude separating Western Accumulation from Pacific
#'   Expansion (default -74).
#' @param surface_only Logical. If TRUE, restricts the analysis to UHW/PMW or
#'   Depth <= 50 m. If FALSE (default), the full water column is used.
#'
#' @return Annual dataframe with:
#'   - n_unique_stations, Missed_West, Missed_East
#'   - Dominance_Index: fraction of samples with fPW >= fPW_cutoff (0.50)
#'   - Front_Lon: easternmost longitude with fPW >= front_detection_thr (0.40);
#'     -999 if no front is detected
#'   - Regime, QC_Flag, Regime_Order (1-6; 0 = insufficient data)
#'
#' @details
#' REGIME CLASSIFICATION (evaluated in this order):
#' 1. Atlantic Dominance      : DI < 0.05
#' 2. Marginal Pacific Inflow : 0.05 <= DI < 0.25
#' 3. Mixed Dominance         : 0.25 <= DI < 0.40
#' 6. Pacific Saturation      : DI >= 0.60 (whatever the front position)
#' 4. Western Accumulation    : 0.40 <= DI < 0.60 AND Front_Lon <= -74°W
#' 5. Pacific Expansion       : 0.40 <= DI < 0.60 AND Front_Lon >  -74°W
#'
#' CHEMICAL DEFINITIONS:
#' - PW-dominated sample : fPW >= 0.50 (PW is the majority end-member)
#' - PW leading edge     : fPW >= 0.40 (diluted PW signal in the frontal zone)
#'
#' SPATIAL COVERAGE:
#' - Western coverage: samples west of -76.0°W
#' - Eastern coverage: samples east of -73.0°W
#'
#' DEPTH/WATER-MASS FILTERING:
#' - Optional (surface_only = TRUE): UHW/PMW or Depth <= 50 m.
#' - Default (FALSE): full water column, as used in the manuscript.
#'
#' @author Nicolas Schiffrine, Takuvik Joint International Laboratory
#' @references Yamamoto-Kawai et al. (2008); Jones et al. (2008)
#' @export
pw_regime_calc <- function(data,
                           min_stations        = 6,
                           fPW_cutoff          = 0.50,
                           front_detection_thr = 0.40,
                           western_limit       = -74,
                           surface_only        = FALSE) {
  
  # ── Thresholds ──────────────────────────────────────────────────────────────
  THRES_COLLAPSE   <- 0.05
  THRES_RETRACTION <- 0.25
  THRES_EXTENSION  <- 0.40
  THRES_SATURATION <- 0.60
  
  # ── Optional surface/halocline filter ───────────────────────────────────────
  if (surface_only) {
    data <- data %>%
      filter(classification %in% c("UHW", "PMW") | Depth <= 50)
  }
  
  data %>%
    filter(!is.na(fPW)) %>%                     # NA-safe for both indices
    group_by(year) %>%
    summarise(
      n_unique_stations = n_distinct(Station_Cluster),
      Missed_West = min(decimalLongitude, na.rm = TRUE) > -76.0,
      Missed_East = max(decimalLongitude, na.rm = TRUE) < -73.0,
      
      # Sample-based dominance: fraction of PW-dominated samples
      Dominance_Index = mean(fPW >= fPW_cutoff),
      
      # Leading edge: easternmost longitude with detectable PW
      Front_Lon = {
        lon_front <- decimalLongitude[fPW >= front_detection_thr &
                                        !is.na(decimalLongitude)]
        if (length(lon_front) > 0) max(lon_front) else -999
      },
      .groups = "drop"
    ) %>%
    mutate(
      Regime = case_when(
        n_unique_stations < min_stations      ~ "Insufficient Data",
        Dominance_Index < THRES_COLLAPSE      ~ "Atlantic Dominance",
        Dominance_Index < THRES_RETRACTION    ~ "Marginal Pacific Inflow",
        Dominance_Index < THRES_EXTENSION     ~ "Mixed Dominance",
        Dominance_Index >= THRES_SATURATION   ~ "Pacific Saturation",   # priority
        Front_Lon <= western_limit            ~ "Western Accumulation",
        TRUE                                  ~ "Pacific Expansion"
      ),
      QC_Flag = case_when(
        Regime == "Insufficient Data" ~ "Exclude",
        Missed_West | Missed_East     ~ "Partial Coverage",
        TRUE                          ~ "Good Coverage"
      ),
      Regime_Order = case_when(
        Regime == "Atlantic Dominance"      ~ 1,
        Regime == "Marginal Pacific Inflow" ~ 2,
        Regime == "Mixed Dominance"         ~ 3,
        Regime == "Western Accumulation"    ~ 4,
        Regime == "Pacific Expansion"       ~ 5,
        Regime == "Pacific Saturation"      ~ 6,
        TRUE                                ~ 0
      )
    )
}

# ── Call with explicit thresholds ─────────────────────────────────────────────
# pw_regimes <- pw_regime_calc(
#   pw_regime_df,
#   min_stations        = 6,
#   fPW_cutoff          = 0.50,
#   front_detection_thr = 0.40,
#   western_limit       = -74,
#   surface_only        = FALSE
# )
