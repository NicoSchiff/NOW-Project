# ══════════════════════════════════════════════════════════════════════════════
# GEOGRAPHIC FILTERING FOR NOW TRANSECT
# ══════════════════════════════════════════════════════════════════════════════

#' Extract observations within North Water Polynya spatial bounds
#'
#' Applies geographic filtering to isolate observations from the NOW transect
#' based on predefined latitude and longitude boundaries.
#'
#' @param raw_data data.frame with geographic coordinates
#' @return data.frame filtered to NOW transect region
extract_transect <- function(data) {
  filtered_data <- data %>%
    filter(
      between(decimalLongitude, NOW_BOUNDS$longitude[1], NOW_BOUNDS$longitude[2]),
      between(decimalLatitude, NOW_BOUNDS$latitude[1], NOW_BOUNDS$latitude[2])
    )
  
  # Report filtering results
  cat("Geographic filtering results:\n")
  cat("  Original observations:", nrow(raw_data), "\n")
  cat("  NOW transect observations:", nrow(filtered_data), "\n")
  cat("  Retention rate:", round(100 * nrow(filtered_data) / nrow(raw_data), 1), "%\n")
  
  return(filtered_data)
}
