# ═══════════════════════════════════════════════════════════════════════════════
# ARCTIC WATER MASS CLASSIFICATION USING σ₀-SPICITY METHODOLOGY
# ═══════════════════════════════════════════════════════════════════════════════
#
# This script implements a σ₀-spicity water mass classification approach
# inspired by Huang et al. (2024) for Arctic oceanographic data analysis.
# The methodology provides more robust discrimination of water masses compared 
# to traditional θ-S diagrams, particularly in high-latitude systems.
#
# Author: Nicolas Schiffrine
# Institution: Takuvik Joint International Laboratory, Université Laval
# Associated manuscript: "The North Water Polynya in Transition: Multi-Decadal 
#                        Changes of Hydrographic and Nutrient Dynamics"
# ═══════════════════════════════════════════════════════════════════════════════
# ═══════════════════════════════════════════════════════════════════════════════
# WATER MASS ENDMEMBER DEFINITIONS
# ═══════════════════════════════════════════════════════════════════════════════

#' Define Arctic water mass endmembers
#' 
#' Water mass endmembers based on literature values and regional characteristics
#' for the North Water (NOW) polynya and Baffin Bay system.
#' 
#' @param reference_lon numeric, reference longitude for spicity calculations
#' @param reference_lat numeric, reference latitude for spicity calculations
#' 
#' @details Endmember properties derived from:
#' - Burger et al. (2023) for Arctic water masses
#' - Regional observations from ArcticNet expeditions (2005-2023)
#' - Temperature and salinity uncertainties represent natural variability
#'
#' @return tibble with water mass properties including σ₀ and spicity
#' @export
create_endmembers <- function(reference_lon = -75, reference_lat = 76) {
  
  # Required packages
  library(tidyverse)
  library(purrr)
  library(gsw)      # For seawater properties calculations
  
  
  # Define water mass endmembers with uncertainty estimates
  endmembers <- tibble(
    watermass = c("PMW", "UHW", "CBAW", "BBPW", "SPMW"),
    Temperature = c(-1.5, -1.3, 0.0, -1.1, 0.7),     # Mean conservative temperature (°C)
    Salinity = c(31.2, 32.6, 34.7, 33.6, 34.3),      # Mean absolute salinity (g/kg)
    temp_sd = c(0.1, 0.1, 0.1, 0.3, 0.6),            # Temperature uncertainty (°C)
    sal_sd = c(0.9, 0.3, 0.1, 0.2, 0.1)              # Salinity uncertainty (g/kg)
  ) %>%
    mutate(
      # Calculate potential density anomaly (σ₀) referenced to surface pressure
      sigma0 = gsw_sigma0(SA = Salinity, CT = Temperature),
      
      # Calculate potential spicity (π₀) - orthogonal coordinate to σ₀
      spicity = sw_pspi(
        S = Salinity, 
        temp = Temperature, 
        temp_unit = "conservative",    # Using conservative temperature
        sal_unit = "SA",              # Using absolute salinity
        longitude = reference_lon,     # Reference longitude for calculation
        latitude = reference_lat,      # Reference latitude for calculation
        pr = 0                        # Surface pressure reference
      )
    )
  
  return(endmembers)
}

# ═══════════════════════════════════════════════════════════════════════════════
# MAIN CLASSIFICATION FUNCTION
# ═══════════════════════════════════════════════════════════════════════════════

#' Classify water masses using σ₀-spicity methodology
#'
#' This function implements a water mass classification approach inspired by
#' Huang et al. (2024), using potential density (σ₀) and potential spicity (π₀)
#' coordinates to provide robust water mass discrimination for Arctic waters.
#'
#' @param data data.frame containing oceanographic observations with 'sigma0' and 'spicity' columns
#' @param endmembers data.frame with water mass endmembers (output from create_endmembers())
#' @param threshold numeric, maximum Euclidean distance (kg/m³) for "pure" water mass assignment (default: 0.6)
#' @param comparable_ratio numeric, ratio threshold for detecting transitional waters (default: 1.3)
#' @param return_long logical, return full diagnostic information (TRUE) or just classification (FALSE)
#'
#' @details
#' Classification logic inspired by Huang et al. (2024):
#' 1. Calculate Euclidean distances to all endmembers in σ₀-π₀ space
#' 2. CONDITION 1: Comparable distances to ≥2 endmembers → "Mixed" water
#' 3. CONDITION 2: All distances > threshold → "Mixed" water (far from all endmembers)
#' 4. Otherwise: Assign to closest endmember if distance < threshold
#'
#' @return data.frame with original data plus classification results
#' @export
classify_watermass <- function(data, endmembers, threshold = 0.6, comparable_ratio = 1.3, return_long = FALSE) {
  
  # Input validation
  required_cols <- c("sigma0", "spicity")
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Filter out rows with missing values in required columns
  valid_data <- data %>%
    filter(!is.na(sigma0), !is.na(spicity))
  
  if (nrow(valid_data) == 0) {
    stop("No valid observations found (all sigma0 or spicity values are NA)")
  }
  
  # Perform classification for each valid observation
  results <- valid_data %>%
    rowwise() %>%
    do({
      point <- c(.$sigma0, .$spicity)  # Current observation in σ₀-π₀ space
      
      # Skip if point contains NA values
      if (any(is.na(point))) {
        if (return_long) {
          return(tibble(
            condition_1_comparable = NA,
            condition_2_large_dist = NA,
            classification = NA_character_,
            min_distance = NA_real_,
            closest_endmember = NA_character_,
            PMW = NA_real_, UHW = NA_real_, CBAW = NA_real_, 
            BBPW = NA_real_, SPMW = NA_real_
          ))
        } else {
          return(tibble(classification = NA_character_))
        }
      }
      
      # Calculate Euclidean distances to all endmembers in σ₀-spicity space
      distances <- endmembers %>%
        mutate(distance = sqrt((sigma0 - point[1])^2 + (spicity - point[2])^2)) %>%
        pull(distance, name = watermass)
      
      # CONDITION 1: Comparable distances to multiple endmembers (transitional waters)
      sorted_dist <- sort(distances)
      condition_1 <- length(sorted_dist) >= 2 && 
        !is.na(sorted_dist[1]) && !is.na(sorted_dist[2]) &&
        sorted_dist[1] > 0 &&
        (sorted_dist[2] / sorted_dist[1]) <= comparable_ratio && 
        length(distances[distances <= sorted_dist[2] * comparable_ratio]) >= 2
      
      # CONDITION 2: All distances exceed threshold (complex mixing)
      condition_2 <- all(distances > threshold, na.rm = TRUE)
      
      # Apply classification logic with explicit handling of edge cases
      min_dist <- min(distances, na.rm = TRUE)
      closest_wm <- names(which.min(distances))
      
      classification <- if (is.na(min_dist) || length(closest_wm) == 0) {
        "Unclassified"
      } else if (condition_1 || condition_2) {
        "Mixed"
      } else if (min_dist < threshold) {
        closest_wm
      } else {
        "Unclassified"
      }
      
      # Return results based on requested detail level
      if (return_long) {
        tibble(
          condition_1_comparable = condition_1,
          condition_2_large_dist = condition_2,
          classification = classification,
          min_distance = min_dist,
          closest_endmember = ifelse(length(closest_wm) > 0, closest_wm, NA_character_),
          PMW = distances["PMW"],
          UHW = distances["UHW"],
          CBAW = distances["CBAW"],
          BBPW = distances["BBPW"],
          SPMW = distances["SPMW"]
        )
      } else {
        tibble(
          classification = classification,
          distance = min_dist)
      }
    }) %>%
    ungroup()
  
  # Bind results back to original data
  return(bind_cols(valid_data, results))
}

# ═══════════════════════════════════════════════════════════════════════════════
# USAGE EXAMPLES (commented out for sourcing)
# ═══════════════════════════════════════════════════════════════════════════════

# # RECOMMENDED WORKFLOW:
# 
# # 1. Create endmembers
# endmembers <- create_endmembers(reference_lon = -75, reference_lat = 76)
# 
# # 2. Perform sensitivity analysis
# ratio_analysis <- test_ratios(your_data, endmembers)
# threshold_analysis <- test_thresholds(your_data, endmembers)
# 
# # 3. Apply classification
# classified_data <- classify_watermass(
#   data = your_data,
#   endmembers = endmembers,
#   threshold = 0.6,        # Optimized value
#   comparable_ratio = 1.3, # Balanced transitional water detection
#   return_long = FALSE     # Simple output
# )
# 
# # 4. Check results
# table(classified_data$classification)