query_wm_stat <- function(data, wm = NULL, station = NULL, var = NULL) {
  
  out <- data
  
  if (!is.null(wm))      out <- dplyr::filter(out, classification   %in% wm)
  if (!is.null(station)) out <- dplyr::filter(out, Station_Cluster  %in% station)
  if (!is.null(var))     out <- dplyr::filter(out, variable         %in% var)
  
  if (nrow(out) == 0L) {
    warning("No rows match the requested combination.", call. = FALSE)
    return(out)
  }
  
  out %>%
    dplyr::arrange(Station_Cluster, classification, variable) %>%
    dplyr::select(
      Station_Cluster, classification, variable,
      dplyr::any_of(c("n","min", "max", "median", "iqr"))
    )
}

query_wm_stat(data = wm_stat, wm = "PMW", var = "Salinity")

#' Range agrégé sur toutes les stations pour une masse d'eau
query_wm_range_global <- function(data, wm = NULL, var = NULL, digits = 2) {
  
  out <- data
  if (!is.null(wm))  out <- dplyr::filter(out, classification %in% wm)
  if (!is.null(var)) out <- dplyr::filter(out, variable       %in% var)
  
  out %>%
    dplyr::group_by(classification, variable) %>%
    dplyr::summarise(
      min = min(median, na.rm = TRUE),
      max = max(median, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(range = sprintf("%.*f-%.*f", digits, min, digits, max))
}

query_wm_range_global(wm_stat, wm = "BBPW")
query_wm_stat(wm_stat, wm = "SPMW") %>% 
  pivot_longer(c(n, min, max, median, iqr),
               names_to = "resultat", values_to = "val") %>%
  pivot_wider(names_from = variable, values_from = val) %>% print(n = Inf, width = Inf)


