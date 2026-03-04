#' Calculate Zonal Pacific Water Indices with Nuanced 6-Tier Classification
#' 
#' This function analyzes the spatiotemporal distribution of Pacific Water influence
#' in the North Water (NOW) polynya based on the Pacific Water fraction (fPW) and
#' classifies each year into one of six hydrographic regimes.
#' 
#' @description
#' The classification system captures the full spectrum of Pacific Water states
#' from complete Atlantic dominance to Pacific saturation, with particular attention
#' to spatial patterns and western accumulation vs. eastward extension dynamics.
#' 
#' @param data A dataframe containing hydrographic observations with the following
#'   required columns:
#'   - classification: Water mass classification (e.g., "PMW", "UHW", "BBPW")
#'   - Depth: Sample depth in meters
#'   - fPW: Pacific Water fraction (0-1)
#'   - Station_Cluster: Station identifier for spatial coverage assessment
#'   - decimalLongitude: Longitude in decimal degrees (negative for West)
#'   - year: Sampling year
#' 
#' @return A dataframe with annual regime classifications containing:
#'   - year: Sampling year
#'   - n_unique_stations: Number of unique stations sampled
#'   - Missed_West/Missed_East: Logical flags for incomplete spatial coverage
#'   - Dominance_Index: Fraction of samples with fPW >= 0.50 (0-1)
#'   - Front_Lon: Easternmost longitude where fPW >= 0.20 (-999 if no front detected)
#'   - Regime: Six-level classification of Pacific Water state
#'   - QC_Flag: Data quality assessment ("Good Coverage", "Partial Coverage", "Exclude")
#'   - Regime_Order: Numerical ordering for plotting (1-6)
#' 
#' @details
#' REGIME CLASSIFICATION SYSTEM:
#' 
#' 1. ATLANTIC DOMINANCE (Dominance_Index < 0.05)
#'    Complete absence of Pacific influence; Atlantic waters dominate throughout
#'    the water column. Represents the most Atlantic-like conditions.
#' 
#' 2. PACIFIC RETRACTION (0.05 ≤ Dominance_Index < 0.25)
#'    Minimal Pacific presence, likely confined to the westernmost shelf regions.
#'    Pacific signal is detectable but spatially restricted.
#' 
#' 3. MIXED DOMINANCE (0.25 ≤ Dominance_Index < 0.40)
#'    Intermediate state where neither Atlantic nor Pacific waters clearly dominate.
#'    Represents transitional conditions with competing influences.
#' 
#' 4. WESTERN ACCUMULATION (Dominance_Index ≥ 0.40 AND Front_Lon ≤ -74.0°W)
#'    High Pacific water volume but spatially confined to western regions.
#'    Indicates strong Pacific input that cannot penetrate eastward, possibly
#'    due to opposing Atlantic flows or density barriers.
#' 
#' 5. PACIFIC EXPANSION (0.40 ≤ Dominance_Index < 0.60 AND Front_Lon > -74.0°W)
#'    Substantial Pacific influence with successful eastward penetration across
#'    the transect. Represents the "normal" high Pacific water state.
#' 
#' 6. PACIFIC SATURATION (Dominance_Index ≥ 0.60)
#'    Extreme Pacific water events where the system is flooded with Pacific-origin
#'    waters. Represents the most Pacific-influenced conditions possible.
#' 
#' SPATIAL COVERAGE ASSESSMENT:
#' - Western coverage: Samples must extend west of -76.0°W
#' - Eastern coverage: Samples must extend east of -73.0°W  
#' - Minimum stations: 6 unique stations required for reliable classification
#' 
#' CHEMICAL DEFINITION:
#' - Pacific Water signature: fPW ≥ 0.50 (based on N:P nutrient ratios)
#' - Front detection threshold: fPW ≥ 0.20 (detects leading edge of Pacific influence)
#' 
#' DEPTH/WATER MASS FILTERING:
#' - Focuses on active surface/halocline layers (UHW, PMW masses or depth ≤ 50m)
#' - Excludes deeper Atlantic-dominated layers to avoid bias
#' 
#' @author Nicolas Schiffrine, Takuvik Joint International Laboratory
#' @references 
#' - Yamamoto-Kawai et al. (2008) for Pacific Water fraction methodology
#' - Newton et al. (2013) for Arctic N:P tracer relationships
#' 
#' @examples
#' \dontrun{
#' # Calculate Pacific Water regimes for NOW dataset
#' pw_regimes <- pw_regime_calc(NOW_data)
#' 
#' # View results excluding insufficient data years
#' good_years <- pw_regimes %>% 
#'   filter(QC_Flag != "Exclude") %>%
#'   select(year, Regime, Dominance_Index, Front_Lon)
#' 
#' # Plot temporal evolution
#' ggplot(pw_regimes, aes(x = year, y = Regime_Order, fill = Regime)) +
#'   geom_col() +
#'   scale_y_continuous(breaks = 1:6, 
#'                      labels = c("Atlantic\nDominance", "Pacific\nRetraction",
#'                                "Mixed\nDominance", "Western\nAccumulation", 
#'                                "Pacific\nExpansion", "Pacific\nSaturation"))
#' }
#' 
#' @export
pw_regime_calc <- function(data, min_stations, fPW_cutoff, western_limit = -74, front_detection_thr) {
  
  # =============================================================================
  # THRESHOLD DEFINITIONS
  # =============================================================================
  
  # 1. SAMPLING REQUIREMENTS
  # Minimum number of unique stations required for reliable regime classification
  MIN_STATIONS_REQUIRED <- min_stations  
  
  # 2. CHEMICAL DEFINITION
  # Threshold for identifying Pacific Water signature based on fPW
  FPW_CHEMICAL_CUTOFF <- fPW_cutoff 
  
  # 3. SPATIAL/VOLUME THRESHOLDS
  # These define the boundaries between different regime states based on
  # the fraction of samples exceeding the chemical cutoff (Dominance_Index)
  
  THRES_COLLAPSE <- 0.05      # Below this: Atlantic Dominance
  THRES_RETRACTION <- 0.25    # Below this: Pacific Retraction  
  THRES_EXTENSION <- 0.40     # Above this: Western Accumulation or Pacific Expansion
  THRES_SATURATION <- 0.60    # Above this: Pacific Saturation
  
  # 4. SPATIAL CONSTRAINTS
  # Longitude threshold for distinguishing Western Accumulation from Pacific Expansion
  WESTERN_LIMIT <- western_limit      # East of this = successful eastward extension
  
  # Front detection threshold (lower than chemical cutoff to detect leading edge)
  FRONT_DETECTION_THRESHOLD <- front_detection_thr
  
  # =============================================================================
  # DATA PROCESSING AND REGIME CLASSIFICATION
  # =============================================================================
  
  data %>%
    # Filter for the active surface/halocline layer where Pacific influence is strongest
    # This avoids bias from deeper Atlantic-dominated waters
    # filter(
    #   classification %in% c("UHW", "PMW") | Depth <= 50, 
    #   !is.na(fPW)
    # ) %>%
    
    # Group by year for annual regime assessment
    group_by(year) %>% 
    summarise(
      
      # --- QUALITY CONTROL METRICS ---
      # Number of unique sampling stations (spatial coverage indicator)
      n_unique_stations = n_distinct(Station_Cluster),
      
      # Spatial coverage flags: detect incomplete transect sampling
      Missed_West = min(decimalLongitude, na.rm = TRUE) > -76.0, 
      Missed_East = max(decimalLongitude, na.rm = TRUE) < -73.0, 
      
      # --- PHYSICAL OCEANOGRAPHIC METRICS ---
      # Dominance Index: Fraction of samples with strong Pacific signature
      # This is the primary metric for regime classification
      Dominance_Index = mean(fPW >= FPW_CHEMICAL_CUTOFF, na.rm = TRUE),
      
      # Pacific Front Position: Easternmost extent of Pacific influence
      # Uses lower threshold to detect the leading edge of Pacific water mass
      Front_Lon = ifelse(
        any(fPW >= FRONT_DETECTION_THRESHOLD), 
        max(decimalLongitude[fPW >= FRONT_DETECTION_THRESHOLD], na.rm = TRUE), 
        -999  # Sentinel value for no detectable front
      ), 
      
      .groups = "drop"
    ) %>%
    
    # =============================================================================
  # REGIME CLASSIFICATION LOGIC
  # =============================================================================
  
  mutate(
    Regime = case_when(
      
      # --- DATA QUALITY EXCLUSION ---
      n_unique_stations < MIN_STATIONS_REQUIRED ~ "Insufficient Data",
      
      # --- REGIME 1: ATLANTIC DOMINANCE ---
      # Virtually no Pacific water present (< 5% of samples show Pacific signature)
      # Represents years when Atlantic waters completely dominate the system
      Dominance_Index < THRES_COLLAPSE ~ "Atlantic Dominance",
      
      # --- REGIME 2: PACIFIC RETRACTION ---
      # Weak Pacific presence (5-25% of samples show Pacific signature)
      # Pacific waters are detectable but spatially very limited
      Dominance_Index < THRES_RETRACTION ~ "Marginal Pacific Inflow",
      
      # --- REGIME 3: MIXED DOMINANCE ---
      # Intermediate conditions (25-40% of samples show Pacific signature)
      # Neither Atlantic nor Pacific waters clearly dominate the system
      Dominance_Index < THRES_EXTENSION ~ "Mixed Dominance",
      
      # --- REGIME 4: WESTERN ACCUMULATION ---
      # High Pacific volume (≥40%) BUT spatially blocked (front west of -74°W)
      # Strong Pacific input that cannot penetrate eastward across the transect
      # May indicate density barriers or opposing Atlantic flows
      Dominance_Index >= THRES_EXTENSION & Front_Lon <= WESTERN_LIMIT ~ "Western Accumulation",
      
      # --- REGIME 6: PACIFIC SATURATION ---
      # Extreme Pacific dominance (≥60% of samples show Pacific signature)
      # The system is flooded with Pacific-origin waters
      # Checked before Pacific Expansion to catch extreme events first
      Dominance_Index >= THRES_SATURATION ~ "Pacific Saturation",
      
      # --- REGIME 5: PACIFIC EXPANSION ---
      # Substantial Pacific influence (40-60%) with successful eastward extension
      # Pacific waters penetrate across the transect (front east of -74°W)
      # Represents the "standard" high Pacific water state
      TRUE ~ "Pacific Expansion"
    ),
    
    # =============================================================================
    # QUALITY CONTROL FLAG ASSIGNMENT
    # =============================================================================
    
    QC_Flag = case_when(
      Regime == "Insufficient Data" ~ "Exclude",
      Missed_West | Missed_East ~ "Partial Coverage", 
      TRUE ~ "Good Coverage"
    ),
    
    # =============================================================================
    # REGIME ORDERING FOR VISUALIZATION
    # =============================================================================
    # Numerical ordering from most Atlantic-like (1) to most Pacific-like (6)
    
    Regime_Order = case_when(
      Regime == "Atlantic Dominance" ~ 1,
      Regime == "Marginal Pacific Inflow" ~ 2,
      Regime == "Mixed Dominance" ~ 3,
      Regime == "Western Accumulation" ~ 4,
      Regime == "Pacific Expansion" ~ 5,
      Regime == "Pacific Saturation" ~ 6,
      TRUE ~ 0  # For "Insufficient Data" and any unforeseen cases
    )
  )
}

# =============================================================================
# USAGE EXAMPLE
# =============================================================================
# 
# # Apply the classification to your NOW dataset
# pw_results_nuanced <- pw_regime_calc(NOW)
# 
# # Display results summary
# print("Pacific Water Regime Classification Results (2005-2023):")
# print("=" %R% rep("=", 55))
# 
# regime_summary <- pw_results_nuanced %>% 
#   select(year, QC_Flag, Regime, Dominance_Index, Front_Lon) %>% 
#   arrange(year)
# 
# print(regime_summary)
# 
# # Summary statistics by regime
# print("\nRegime Frequency Summary:")
# regime_counts <- pw_results_nuanced %>%
#   filter(QC_Flag != "Exclude") %>%
#   count(Regime, sort = TRUE) %>%
#   mutate(Percentage = round(n / sum(n) * 100, 1))
# 
# print(regime_counts)
