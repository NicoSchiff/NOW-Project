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
rm(list = ls())  # Clear workspace for reproducible analysis

# Set working directory for nutrient database access
setwd(dir = "/Users/nicolas/Library/CloudStorage/OneDrive-UniversitéLaval/")

# ══════════════════════════════════════════════════════════════════════════════
# CONFIGURATION AND PARAMETERS
# ══════════════════════════════════════════════════════════════════════════════

# Define biogeochemical parameters of interest for NOW analysis
parameter <- c("Temperature", "Salinity", "Nitrate", "Nitrite", "Ammonium", "Oxconc",
               "Phosphate", "Silicate")

# Spatial boundaries for North Water Polynya region (Baffin Bay/Nares Strait)
# Based on oceanographic literature and regional bathymetry
NOW_BOUNDS <- list(
  longitude = c(-80, -65),  # Western to Eastern extent (degrees West)
  latitude = c(76.2, 76.5)  # Southern to Northern extent (degrees North)
)

# ══════════════════════════════════════════════════════════════════════════════
# EXTERNAL FUNCTION LOADING
# ══════════════════════════════════════════════════════════════════════════════

# Source specialized functions for oceanographic analysis
source("DataBase_Nut/20260209/20260209_Package.R")     # σ₀-spicity water mass classification
source("DataBase_Nut/20260209/20260209_DataLoading.R")     # σ₀-spicity water mass classification
source("DataBase_Nut/20260209/20260209_Extract_Transect.R")     # σ₀-spicity water mass classification
source("DataBase_Nut/20260209/20260209_ClusterStation.R")     # σ₀-spicity water mass classification
source("DataBase_Nut/20260209/20260209_sw_pspi.r")               # Seawater spicity calculations
source("DataBase_Nut/20260209/20260209_WaterMass_Classif.R")     # σ₀-spicity water mass classification
source("DataBase_Nut/20260209/20260209_PacificWater.R")     # σ₀-spicity water mass classification

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
  n_clusters = 15,                    # Target number based on transect length
  min_stations_per_cluster = 7        # Ensure sufficient temporal coverage per cluster
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
            by = c("Station", "decimalLongitude", "decimalLatitude", "year", "month", "day")) %>%
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
)

cat("Water mass endmembers created:\n")
print(wm_ref %>% select(watermass, Temperature, Salinity, sigma0, spicity))

# PARAMETER OPTIMIZATION ANALYSIS
# Note: This requires NOW.WaterMass to be defined first
# The following analysis should be run after calculating oceanographic properties

now_hydro <- 
  now_geo %>% 
  mutate(
    # Calculate potential density anomaly (σ₀) referenced to surface pressure
    sigma0 = gsw_sigma0(SA = Salinity, CT = Temperature),
    
    # Calculate potential spicity (π₀) - orthogonal coordinate to σ₀
    # Spicity quantifies water mass mixing independent of density
    spicity = sw_pspi(
      S = Salinity, 
      temp = Temperature, 
      temp_unit = "conservative",    # Using conservative temperature
      sal_unit = "SA",              # Using absolute salinity
      longitude = ref_coords$longitude,     # Reference longitude for calculation
      latitude = ref_coords$latitude,      # Reference latitude for calculation
      pr = 0                        # Surface pressure reference
    )
  )

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
  wm_ref,
  threshold = 0.6,        # Optimized based on Huang et al. (2024)
  comparable_ratio = 1.3, # Balanced transitional water detection
  return_long = F     # Simple output for routine analysis
)

# # Generate classification summary statistics
# classification_summary <- now_wm %>%
#   count(classification) %>%
#   arrange(desc(n)) %>%
#   mutate(percentage = round(100 * n / sum(n), 1))
# 
# cat("=== WATER MASS CLASSIFICATION RESULTS ===\n")
# print(classification_summary)

# ══════════════════════════════════════════════════════════════════════════════
# DERIVED VARIABLES
# ══════════════════════════════════════════════════════════════════════════════
NOW <-
  now_wm %>%
  mutate(
    # >>> “N*” after Gruber & Sarmiento.
    Nstar = 0.87 * (Nitrate - 16 * Phosphate + 2.9),
    
    # >>> Total inorganic nitrogen: measured NO3 + NO2 and, if present, NH4.
    #     `exists()` protects the pipeline in case “Ammonium” is not in the file.
    TIN = Nitrate + 
      replace_na(Nitrite, 0) + 
      replace_na(Ammonium, 0),
    
    # >>> “Arctic N–P” tracer (ANP) following Newton et al. 2013.
    ANP = ANP(Phosphate = Phosphate, Nitrogen = TIN),
    
    # >>> Fraction of Pacific Water (%). Function returns 0-1, here multiplied by 100.
    fPW = fpw(
      Phosphate = Phosphate, Nitrogen = TIN,
      atlantic_ref = "Jones1998",
      pacific_ref = "Yamamoto.Kawai2008"
    ),
    
    # >>> Semiconservative tracer combining nitrate, Phosphate and dissolved oxygen.
    NO = (9 * Nitrate) + Oxconc,
    PO = (135 * Phosphate) + +Oxconc,
    NO_PO = NO / PO,
    
    # >>> Salinity-dependent POs* tracer for Arctic waters
    POs_star = POs_star(phosphate = Phosphate, oxygen = Oxconc, salinity = Salinity),
    classification = factor(classification, c("PMW","UHW", "BBPW", "SPMW", "CBAW", "Mixed")),
    
    # >>> Station (cluster) order west → east.
    Station_Cluster = factor(Station_Cluster, c("100", "101", "103", "105", "107", "108", "110", "111", "113", "115", "116"))
    
  ) %>%
  select(
    Cruise:Salinity,sigma0:spicity,
    Nitrate:Silicate,TIN,Nstar,NO:POs_star,
    classification,distance,ANP,fPW
  ) %>% 
  filter(!year %in% c(1997, 1998, 1999, 2007, 2008))

