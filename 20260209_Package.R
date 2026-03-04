# ══════════════════════════════════════════════════════════════════════════════
# REQUIRED PACKAGES FOR NOW TRANSECT ANALYSIS
# ══════════════════════════════════════════════════════════════════════════════

required_packages <- c(
  "broom", "colorspace", "factoextra", "fields", "forcats", 
  "GGally", "ggcorrplot", "ggh4x", "ggnewscale", "ggoce", 
  "ggOceanMaps", "ggplot2", "ggpmisc", "ggpubr", "ggrepel", 
  "ggspatial", "grid", "gridExtra", "gsw", "lubridate", 
  "MBA", "metR", "pals", "RColorBrewer", "rstatix", 
  "see", "tibble", "tidyverse", "trend"
)


#' Load all required packages for NOW analysis
#' @param packages vector of package names
#' @param install_missing logical, install if missing
load_now_packages <- function(packages = required_packages, install_missing = TRUE) {
  
  for(pkg in packages) {
    if(!require(pkg, character.only = TRUE, quietly = TRUE)) {
      if(install_missing) {
        install.packages(pkg, dependencies = TRUE)
        library(pkg, character.only = TRUE)
      } else {
        warning("Package ", pkg, " not found and install_missing = FALSE")
      }
    }
  }
  
  cat("Loaded", length(packages), "packages for NOW analysis ✅\n")
}

# Charger tout
load_now_packages(required_packages)
