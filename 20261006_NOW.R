# ══════════════════════════════════════════════════════════════════════════════
# NORTH WATER POLYNYA (NOW) TRANSECT DATA PROCESSING AND STATION CLUSTERING
# ══════════════════════════════════════════════════════════════════════════════
#
# Comprehensive data processing pipeline for the North Water Polynya (NOW) transect
# analysis covering multi-decadal hydrographic and biogeochemical observations
# from ArcticNet expeditions (2005-2023).
#
# This script performs:
# 1. Multi-year nutrient database loading and standardization
# 2. Geographic filtering for NOW transect extraction
# 3. Hierarchical station clustering for spatial consistency
# 4. Water mass classification using σ₀-spicity methodology
# 5. Validation and visualization of results
#
# Author: Nicolas Schiffrine
# Institution: Takuvik Joint International Laboratory, Université Laval
# Associated manuscript: "The North Water Polynya in Transition: Multi-Decadal
#                        Changes of Hydrographic and Nutrient Dynamics"
# ══════════════════════════════════════════════════════════════════════════════

# Environment initialization -----------------------------------------------
rm(list = ls()) # Clear workspace for reproducible analysis

# Set working directory for nutrient database access
setwd(dir = "/Users/nicolas/Library/CloudStorage/OneDrive-UniversitéLaval/")

# ══════════════════════════════════════════════════════════════════════════════
# CONFIGURATION AND PARAMETERS
# ══════════════════════════════════════════════════════════════════════════════

# Define biogeochemical parameters of interest for NOW analysis
parameter <- c(
  "Temperature", "Salinity", "Nitrate", "Nitrite", "Ammonium", "Oxconc",
  "Phosphate", "Silicate"
)

nut_vars <- c("Nitrate", "Nitrite", "Ammonium", "Oxconc", "Phosphate", "Silicate")

# Spatial boundaries for North Water Polynya region (Baffin Bay/Nares Strait)
# Based on oceanographic literature and regional bathymetry
NOW_BOUNDS <- list(
  longitude = c(-80, -65), # Western to Eastern extent (degrees West)
  latitude = c(76.2, 76.5) # Southern to Northern extent (degrees North)
)

# ══════════════════════════════════════════════════════════════════════════════
# EXTERNAL FUNCTION LOADING
# ══════════════════════════════════════════════════════════════════════════════

# Source specialized functions for oceanographic analysis
source("DataBase_Nut/20260209/20260209_Package.R") # σ₀-spicity water mass classification
source("DataBase_Nut/20260209/20260209_DataLoading.R") # σ₀-spicity water mass classification
source("DataBase_Nut/20260209/20260209_Extract_Transect.R") # σ₀-spicity water mass classification
source("DataBase_Nut/20260209/20260209_ClusterStation.R") # σ₀-spicity water mass classification
source("DataBase_Nut/20260209/20260209_sw_pspi.r") # Seawater spicity calculations
source("DataBase_Nut/20260209/20260209_WaterMass_Classif.R") # σ₀-spicity water mass classification
source("DataBase_Nut/20260209/20260209_PacificWater.R") # σ₀-spicity water mass classification
source("DataBase_Nut/20260209/20260209_ratio_calc.R")

# ══════════════════════════════════════════════════════════════════════════════
# DATA LOADING AND PREPROCESSING
# ══════════════════════════════════════════════════════════════════════════════
# Execute data loading pipeline
raw_data <- load_data()

# ══════════════════════════════════════════════════════════════════════════════
# GEOGRAPHIC FILTERING FOR NOW TRANSECT
# ══════════════════════════════════════════════════════════════════════════════
# Apply spatial filtering to isolate NOW transect data
raw_now <- extract_transect(data = raw_data)

# ══════════════════════════════════════════════════════════════════════════════
# STATION CLUSTERING FOR SPATIAL CONSISTENCY
# ══════════════════════════════════════════════════════════════════════════════
# Apply hierarchical clustering with optimized parameters
now_clusters <- cluster_stations(
  data = raw_now,
  n_clusters = 15, # Target number based on transect length
  min_stations_per_cluster = 7 # Ensure sufficient temporal coverage per cluster
)

# ══════════════════════════════════════════════════════════════════════════════
# CLUSTER NAMING WITH SPATIAL CONFLICT RESOLUTION
# ══════════════════════════════════════════════════════════════════════════════
# Apply improved cluster naming with spatial conflict resolution
cat("=== CLUSTER NAMING WITH CONFLICT RESOLUTION ===\n")
now_renamed <- rename_clusters(now_clusters)

# # Merge clustered station assignments with full NOW dataset
now_geo <- raw_now %>%
  left_join(now_renamed,
    by = c("Station", "decimalLongitude", "decimalLatitude", "year", "month", "day")
  ) %>%
  # Reorganize columns for analytical clarity
  select(Cruise, Station, Station_Cluster, everything(), -Cluster) %>%
  # Retain only observations assigned to valid clusters
  filter(!is.na(Station_Cluster))

# # Generate final dataset summary
# cat("=== FINAL NOW TRANSECT DATASET SUMMARY ===\n")
# cat("Total observations after clustering:", nrow(now_geo), "\n")
# cat("Valid station clusters:", length(unique(now_geo$Station_Cluster)), "\n")
# cat("Temporal coverage:", paste(range(now_geo$year, na.rm = TRUE), collapse = " - "), "\n")

# ══════════════════════════════════════════════════════════════════════════════
# WATER MASS CLASSIFICATION ANALYSIS
# ══════════════════════════════════════════════════════════════════════════════

cat("=== WATER MASS CLASSIFICATION USING σ₀-SPICITY METHODOLOGY ===\n")

# Calculate reference coordinates once at the beginning
ref_coords <- list(
  longitude = median(now_geo$decimalLongitude, na.rm = TRUE),
  latitude = median(now_geo$decimalLatitude, na.rm = TRUE)
)

# Create water mass endmembers for classification
wm_ref <- create_endmembers(
  reference_lon = ref_coords$longitude,
  reference_lat = ref_coords$latitude
) %>% 
  filter(!watermass %in% "CBAW")

cat("Water mass endmembers created:\n")
print(wm_ref %>% select(watermass, CT, SA, sigma0, spicity))

# PARAMETER OPTIMIZATION ANALYSIS
# Note: This requires NOW.WaterMass to be defined first
# The following analysis should be run after calculating oceanographic properties

now_hydro <-
  now_geo %>%
  mutate(
    # 1. Sea pressure (dbar) from depth (m); z is negative downward in TEOS-10
    Pressure = gsw_p_from_z(z = -Depth, latitude = decimalLatitude),
    
    # 2. Practical Salinity (SP, PSS-78, unitless) -> Absolute Salinity (SA, g kg-1)
    #    Station coordinates are used because the Absolute Salinity Anomaly (δSA) varies spatially
    SA = gsw_SA_from_SP(
      SP = Salinity, p = Pressure,
      longitude = decimalLongitude, latitude = decimalLatitude
    ),
    
    # 3. In-situ temperature (t, ITS-90, °C) -> Conservative Temperature (CT, °C)
    CT = gsw_CT_from_t(SA = SA, t = Temperature, p = Pressure),
    
    # 4. Potential temperature referenced to 0 dbar (θ, °C)
    #    Needed for the O2 solubility function (step 7)
    pt = gsw_pt_from_CT(SA = SA, CT = CT),
    
    # 5. Potential density anomaly referenced to 0 dbar (σ0 = ρ(SA, CT, 0) - 1000, kg m-3)
    sigma0 = gsw_sigma0(SA = SA, CT = CT),
    
    # 6. In-situ density (kg m-3 / 1000 = kg L-1); includes pressure compression
    rho = gsw_rho(SA = SA, CT = CT, p = Pressure) / 1000,
    
    # 7. Potential spicity referenced to 0 dbar (π0, kg m-3)
    #    Computed from the same SA and CT as sigma0, so it is orthogonal to σ0 isopycnals
    spicity = gsw_pspi(SA = SA, CT = CT, pr = 0),
    
    # 8. O2 solubility at 1 atm (µmol kg-1; Garcia & Gordon, 1992)
    #    Requires Practical Salinity and potential temperature
    # oxsat = gsw_O2sol_SP_pt(SP = Salinity, pt = pt),
    
    # 9. Apparent Oxygen Utilisation (AOU = O2,sat - O2,obs, µmol kg-1)
    #    Oxconc must already be in µmol kg-1
    # aou = oxsat - Oxconc
  )

# Remove samples without valid T/S (cannot be classified)
n_before <- nrow(now_hydro)

now_hydro <- now_hydro %>%
  filter(!is.na(Salinity), !is.na(Temperature),
         !is.na(sigma0), !is.na(spicity))

cat("Removed", n_before - nrow(now_hydro), "samples without valid T/S;",
    nrow(now_hydro), "samples retained.\n")

# Sanity check: no NA left in the variables used for classification
stopifnot(!anyNA(now_hydro[, c("SA", "CT", "sigma0", "spicity")]))

# ──────────────────────────────────────────────────────────────────────────────
# RUN ANALYSIS
# ──────────────────────────────────────────────────────────────────────────────

# source("DataBase_Nut/20260209/20260209_WaterMass_Sensitivity.R")     # σ₀-spicity water mass classification
# 
# wm_sensitivity_data <- now_hydro %>%
#   select(contains("sigma"), contains("spic"), contains("density"))
# 
# sensitivity_results <- test_parameter_grid(data = wm_sensitivity_data, endmembers = wm_ref)
# performance_metrics <- calculate_performance_metrics(grid_results = sensitivity_results)

# Apply optimized classification to NOW dataset
now_wm <- classify_watermass(
  now_hydro,
  wm_ref %>% filter(!watermass %in% "CBAW"),
  threshold = 0.6, # Optimized based on Huang et al. (2024)
  comparable_ratio = 1.3, # Balanced transitional water detection
  return_long = F # Simple output for routine analysis
)

# # Generate classification summary statistics
classification_summary <- now_wm %>%
  count(classification) %>%
  arrange(desc(n)) %>%
  mutate(percentage = round(100 * n / sum(n), 1))

cat("=== WATER MASS CLASSIFICATION RESULTS ===\n")
print(classification_summary)

# ══════════════════════════════════════════════════════════════════════════════
# DERIVED VARIABLES
# ══════════════════════════════════════════════════════════════════════════════
NOW <-
  now_wm %>%
  mutate(
      across(any_of(nut_vars), ~ .x, .names = "{.col}_L"),
      across(any_of(nut_vars), ~ .x  / rho),
      
    # >>> “N*” after Gruber & Sarmiento.
    Nstar = 0.87 * (Nitrate - 16 * Phosphate + 2.9),

    # >>> Total inorganic nitrogen: measured NO3 + NO2 and, if present, NH4.
    #     `exists()` protects the pipeline in case “Ammonium” is not in the file.
    TIN = Nitrate +
      replace_na(Nitrite, 0) +
      replace_na(Ammonium, 0),

        # >>> Nutrient ratios
    NP = ratio_calc(Nitrate, Phosphate, LOD["Nitrate"], LOD["Phosphate"]),
    SiN = ratio_calc(Silicate, Nitrate, LOD["Silicate"], LOD["Nitrate"]),
    SiP = ratio_calc(Silicate, Phosphate, LOD["Silicate"], LOD["Phosphate"]),

    # >>> “Arctic N–P” tracer (ANP) following Newton et al. 2013.
    ANP = ANP(Phosphate = Phosphate, Nitrogen = TIN),

    # >>> Fraction of Pacific Water (%). Function returns 0-1, here multiplied by 100.
    fPW = fpw(
      Phosphate = Phosphate, Nitrogen = TIN,
      atlantic_ref = "Jones2008",
      pacific_ref = "Yamamoto.Kawai2008"
    ),

    # >>> Semiconservative tracer combining nitrate, Phosphate and dissolved oxygen.
    NO = (9 * Nitrate) + Oxconc,
    PO = (135 * Phosphate) + +Oxconc,
    NO_PO = NO / PO,

    # >>> Salinity-dependent POs* tracer for Arctic waters
    POs_star = POs_star(phosphate = Phosphate, oxygen = Oxconc, salinity = Salinity),
    classification = factor(classification, c("PMW", "UHW", "BBPW", "SPMW", #"CBAW", 
                                              "Mixed")),

    # >>> Station (cluster) order west → east.
    Station_Cluster = factor(Station_Cluster, c("100", "101", "103", "105", "107", "108", "110", "111", "113", "115", "116"))
  ) %>%
  # select(
  #   Cruise:Salinity, CT, SA, sigma0:spicity,
  #   Nitrate:Silicate, TIN,
  #   NP, SiN, SiP,
  #   Nstar, NO:POs_star,
  #   classification, distance, ANP, fPW
  # ) %>%
  filter(!year %in% c(1997, 1998, 1999, 2007, 2008))

