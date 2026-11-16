# ══════════════════════════════════════════════════════════════════════════════
# STATION CLUSTERING FOR SPATIAL CONSISTENCY
# ══════════════════════════════════════════════════════════════════════════════

#' Cluster stations using hierarchical approach for spatial consistency
#'
#' Groups geographically proximate stations to create a consistent spatial framework
#' across multi-year expeditions. This addresses station coordinate variability
#' due to navigation, ice conditions, and operational constraints.
#'
#' @param data data.frame with station coordinates
#' @param n_clusters integer, target number of clusters for hierarchical cutting
#' @param min_stations_per_cluster integer, minimum observations required per cluster
#' @return data.frame with cluster assignments
#'
#' @details
#' Methodology:
#' 1. Extract unique coordinate combinations for clustering analysis
#' 2. Calculate Euclidean distance matrix between all station pairs  
#' 3. Apply Ward's hierarchical clustering to minimize within-cluster variance
#' 4. Cut dendrogram to obtain specified number of clusters
#' 5. Filter clusters with insufficient temporal coverage
cluster_stations <- function(data, n_clusters = 10, min_stations_per_cluster = 5) {
  
  # Extract unique coordinate combinations for clustering analysis
  coordinates <- data %>%
    distinct(Station, decimalLongitude, decimalLatitude, year, month, day)
  
  cat("Station clustering analysis:\n")
  cat("  Unique coordinates for clustering:", nrow(coordinates), "\n")
  
  # Prepare coordinate matrix for distance calculations
  coord_matrix <- coordinates[, c("decimalLatitude", "decimalLongitude")]
  
  # Calculate Euclidean distance matrix between all station pairs
  dist_matrix <- dist(coord_matrix)
  
  # Apply Ward's hierarchical clustering to minimize within-cluster variance
  # Ward's method is optimal for creating compact, spherical clusters
  hc_result <- hclust(dist_matrix, method = "ward.D2")
  
  # Cut dendrogram to obtain specified number of clusters
  cluster_labels <- cutree(hc_result, k = n_clusters)
  
  # Assign cluster identifiers to coordinate data
  coordinates$Cluster <- cluster_labels
  
  # Calculate cluster sizes for threshold filtering
  cluster_counts <- coordinates %>%
    count(Cluster) %>%
    arrange(desc(n))
  
  cat("  Initial clusters created:", n_clusters, "\n")
  cat("  Cluster size distribution:\n")
  print(cluster_counts)
  
  # Identify clusters with insufficient temporal coverage
  small_clusters <- cluster_counts %>%
    filter(n < min_stations_per_cluster) %>%
    pull(Cluster)
  
  # Exclude undersampled clusters from analysis
  coordinates$Cluster[coordinates$Cluster %in% small_clusters] <- NA
  
  # Report clustering results
  cat("  Clusters removed (< ", min_stations_per_cluster, " stations):", length(small_clusters), "\n")
  cat("  Final number of valid clusters:", 
      length(unique(coordinates$Cluster[!is.na(coordinates$Cluster)])), "\n")
  
  return(coordinates)
}

# ══════════════════════════════════════════════════════════════════════════════
# CLUSTER NAMING WITH SPATIAL CONFLICT RESOLUTION
# ══════════════════════════════════════════════════════════════════════════════

#' Assign meaningful names to clusters with spatial conflict resolution
#'
#' Creates interpretable cluster names based on dominant station identity
#' while resolving cases where multiple clusters share the same station name
#' using geographic position-based suffixes.
#'
#' @param cluster_data data.frame with cluster assignments
#' @return data.frame with resolved cluster names
#'
#' @details
#' Conflict resolution strategy:
#' - Single cluster per station: Use station name directly
#' - Multiple clusters per station: Apply W/E/numeric suffixes based on longitude ranking
#' - This maintains oceanographic interpretability while ensuring uniqueness
rename_clusters <- function(cluster_data) {
  
  # Calculate cluster statistics and dominant station for each cluster
  cluster_stats <- cluster_data %>%
    filter(!is.na(Cluster)) %>%
    group_by(Cluster) %>%
    summarise(
      # Identify most frequently sampled station within cluster
      most_frequent_station = names(sort(table(Station), decreasing = TRUE))[1],
      station_frequency = max(table(Station)),
      total_obs = n(),
      
      # Calculate geographic centroid for spatial referencing
      median_lon = median(decimalLongitude),
      median_lat = median(decimalLatitude),
      .groups = "drop"
    )
  
  # Detect naming conflicts where multiple clusters share dominant station name
  station_conflicts <- cluster_stats %>%
    group_by(most_frequent_station) %>%
    summarise(n_clusters = n(), .groups = "drop") %>%
    filter(n_clusters > 1) %>%
    pull(most_frequent_station)
  
  # Resolve conflicts using geographic position-based suffixes
  cluster_names <- cluster_stats %>%
    mutate(
      Station_Cluster = case_when(
        most_frequent_station %in% station_conflicts ~ {
          # Group clusters by station name and rank by longitude
          cluster_order <- ave(median_lon, most_frequent_station, FUN = rank)
          max_clusters <- ave(median_lon, most_frequent_station, FUN = length)
          
          case_when(
            cluster_order == 1 ~ paste0(most_frequent_station, "W"),          # Westernmost
            cluster_order == max_clusters ~ paste0(most_frequent_station, "E"), # Easternmost  
            TRUE ~ paste0(most_frequent_station, "_", cluster_order)          # Intermediate
          )
        },
        TRUE ~ most_frequent_station  # No conflict: use station name directly
      )
    ) %>%
    select(Cluster, Station_Cluster)
  
  # Join cluster names back to original coordinate data
  cluster_data_named <- cluster_data %>%
    left_join(cluster_names, by = "Cluster") %>%
    mutate(Station_Cluster = ifelse(is.na(Cluster), NA, Station_Cluster))
  
  # Report conflict resolution results
  if(length(station_conflicts) > 0) {
    cat("Station naming conflicts detected and resolved for:", 
        paste(station_conflicts, collapse = ", "), "\n")
    conflict_resolution <- cluster_names %>% 
      filter(grepl("W|E|_", Station_Cluster)) %>%
      arrange(Station_Cluster)
    print(conflict_resolution)
  }
  
  return(cluster_data_named)
}


# ══════════════════════════════════════════════════════════════════════════════
# CLUSTER VALIDATION AND QUALITY CONTROL
# ══════════════════════════════════════════════════════════════════════════════

# Validate cluster geographic coherence and generate summary statistics
# cluster_validation <- now_renamed %>%
#   filter(!is.na(Station_Cluster)) %>%
#   group_by(Station_Cluster) %>%
#   summarise(
#     n_obs = n(),                                                    # Total observations per cluster
#     lon_range = paste(round(range(decimalLongitude), 2), collapse = " to "),  # Longitude extent
#     lat_range = paste(round(range(decimalLatitude), 3), collapse = " to "),   # Latitude extent
#     lon_span = diff(range(decimalLongitude)),                      # Cluster geographic span
#     lat_span = diff(range(decimalLatitude)),
#     .groups = "drop"
#   ) %>%
#   arrange(Station_Cluster)
# 
# cat("=== CLUSTER VALIDATION SUMMARY ===\n")
# print(cluster_validation)
# 
# # Quality control check: Flag clusters with excessive geographic spread
# large_clusters <- cluster_validation %>%
#   filter(lon_span > 0.5 | lat_span > 0.1)  # Thresholds based on typical station spacing
# 
# if(nrow(large_clusters) > 0) {
#   cat("\nWarning: Clusters with large geographic spread detected:\n")
#   print(large_clusters)
# }
# 
 
# # Summary by station cluster
# station_summary <- now_geo %>%
#   group_by(Station_Cluster) %>%
#   summarise(
#     observations = n(),
#     year_range = paste(range(year, na.rm = TRUE), collapse = "-"),
#     depth_range = paste(round(range(Depth, na.rm = TRUE)), collapse = "-"),
#     .groups = "drop"
#   ) %>%
#   arrange(Station_Cluster)
# 
# print("Observations per station cluster:")
# print(station_summary)
