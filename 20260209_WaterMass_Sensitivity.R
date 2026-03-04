
# ═══════════════════════════════════════════════════════════════════════════════
# PARAMETER SENSITIVITY ANALYSIS FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════════════

  #' Test sensitivity of classification to comparable ratio parameter
#'
#' @param data data.frame with sigma0 and spicity columns
#' @param endmembers data.frame with water mass endmembers
#' @param ratios numeric vector of comparable ratio values to test
#' @param threshold numeric, fixed threshold for analysis
#' @return data.frame with classification counts for each ratio
#' @export
test_ratios <- function(data, endmembers, ratios = c(1.0, 1.1, 1.2, 1.3, 1.5, 2.0), threshold = 0.6) { 
  
  results <- map_dfr(ratios, ~{
    classified <- classify_watermass(data, endmembers, 
                                     threshold = threshold, comparable_ratio = .x)
    classified %>%
      count(classification) %>%
      mutate(ratio = .x)
  })
  
  return(results %>% 
           pivot_wider(names_from = classification, values_from = n, values_fill = 0))
}

#' Test sensitivity of classification to distance threshold parameter  
#'
#' @param data data.frame with sigma0 and spicity columns
#' @param endmembers data.frame with water mass endmembers
#' @param thresholds numeric vector of threshold values to test
#' @param comparable_ratio numeric, fixed comparable ratio for analysis
#' @return data.frame with classification counts for each threshold
#' @export
test_thresholds <- function(data, endmembers, 
                            thresholds = c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0), 
                            comparable_ratio = 1.3) {
  
  results <- map_dfr(thresholds, ~{
    classified <- classify_watermass(data, endmembers, 
                                     threshold = .x, comparable_ratio = comparable_ratio)
    classified %>%
      count(classification) %>%
      mutate(threshold = .x)
  })
  
  return(results %>% 
           pivot_wider(names_from = classification, values_from = n, values_fill = 0))
}

# ══════════════════════════════════════════════════════════════════════════════
# PARAMETER SENSITIVITY WITH ROBUST NA HANDLING
# ══════════════════════════════════════════════════════════════════════════════

library(tidyverse)
library(viridis)
library(corrplot)

#' Complete 2D parameter grid search with robust column handling
test_parameter_grid <- function(data, endmembers,
                                threshold_seq = seq(0.3, 1.0, 0.1),
                                ratio_seq = seq(1.0, 2.0, 0.2)) {
  
  cat("Testing", length(threshold_seq), "×", length(ratio_seq), "parameter combinations\n")
  
  # Define ALL possible water mass columns
  all_water_masses <- c("BBPW", "CBAW", "Mixed", "PMW", "SPMW", "UHW")
  
  param_grid <- expand_grid(
    threshold = threshold_seq,
    ratio = ratio_seq
  )
  
  results <- map2_dfr(param_grid$threshold, param_grid$ratio, ~{
    
    cat("Testing threshold =", .x, ", ratio =", .y, "\n")
    
    # Run classification
    classified <- classify_watermass(
      data = data,
      endmembers = endmembers, 
      threshold = .x,
      comparable_ratio = .y
    )
    
    # Count water masses with EXPLICIT missing columns handling
    counts <- classified %>%
      count(classification) %>%
      # FORCE all water mass columns to exist with 0 as default
      complete(classification = all_water_masses, fill = list(n = 0)) %>%
      pivot_wider(names_from = classification, values_from = n, values_fill = 0)
    
    # Ensure ALL columns exist (in case complete() missed some)
    for(wm in all_water_masses) {
      if(!wm %in% names(counts)) {
        counts[[wm]] <- 0
      }
    }
    
    # Add parameters and reorder columns
    counts %>% 
      mutate(threshold = .x, ratio = .y, .before = 1) %>%
      select(threshold, ratio, all_of(all_water_masses))
  })
  
  return(results)
}

# Calculate performance metrics with robust NA handling
calculate_performance_metrics <- function(grid_results) {
  
  grid_results %>%
    # Replace any remaining NAs with 0 (shouldn't happen but safety first)
    mutate(across(BBPW:UHW, ~replace_na(.x, 0))) %>%
    
    mutate(
      # Total samples (now guaranteed to work)
      total = BBPW + CBAW + Mixed + PMW + SPMW + UHW,
      
      # Key metrics
      mixed_percent = Mixed / total * 100,
      classification_rate = 100,  # Always 100% since we classify everything
      
      # Shannon diversity index
      diversity = {
        props <- select(., BBPW:UHW) / total
        -rowSums(props * log(props + 1e-10), na.rm = TRUE)  # Add small constant to avoid log(0)
      },
      
      # Balance of pure water masses (CV = coefficient of variation)
      pure_masses = select(., BBPW, PMW, SPMW, UHW),  # Exclude Mixed and rare CBAW
      balance_cv = {
        pure_counts <- select(., BBPW, PMW, SPMW, UHW)
        apply(pure_counts, 1, function(x) sd(x) / mean(x))
      },
      
      # Rare water mass detection (CBAW presence)
      rare_detection = CBAW > 0,
      
      # Classification efficiency (high diversity + balanced + some mixed)
      efficiency = diversity * (1 / (1 + balance_cv)) * (1 - abs(mixed_percent - 15) / 15)
    ) %>%
    
    # Handle any infinite or NaN values
    mutate(
      balance_cv = replace(balance_cv, !is.finite(balance_cv), 0),
      efficiency = replace(efficiency, !is.finite(efficiency), 0)
    )
}

# # Create performance indices
# performance_indices <- grid_metrics %>%
#   mutate(
#     # Normalize metrics to 0-1 scale
#     mixed_score = 1 - abs(mixed_percent - 15) / 15,  # Target: 15% mixed
#     diversity_score = (diversity - min(diversity)) / (max(diversity) - min(diversity)),
#     balance_score = 1 / (1 + balance_cv),  # Lower CV is better
#     rare_score = as.numeric(rare_detection),  # Binary: 0 or 1
#     
#     # Composite indices
#     ecological_index = (diversity_score + balance_score) / 2,
#     technical_index = (mixed_score + rare_score) / 2,
#     performance_index = (ecological_index + technical_index) / 2
#   ) %>%
#   # Clip scores to valid range [0, 1]
#   mutate(across(ends_with("_score"), ~pmax(0, pmin(1, .x))))

