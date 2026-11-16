# ══════════════════════════════════════════════════════════════════════════════
# TEMPORAL TREND VISUALIZATION FUNCTION FOR ARCTIC OCEANOGRAPHIC DATA
# ══════════════════════════════════════════════════════════════════════════════
# 
# Purpose: Generate publication-quality plots showing temporal trends in 
#          oceanographic variables across Arctic water masses
# 
# Author: [Your Name]
# Date: [Current Date]
# 
# Dependencies: tidyverse, ggplot2
# 
# Input data structure expected:
# - Station_Cluster: Station identifier
# - classification: Water mass classification
# - variable: Oceanographic variable name
# - Value: Measured values
# - dec_year: Decimal year for temporal analysis
# - sen_slope: Sen's slope estimator
# - p_value_mk: Mann-Kendall p-value
# - p_value_signif: Statistical significance indicator
# ══════════════════════════════════════════════════════════════════════════════

plot_trend <- function(data, station, water_mass, variable) {
  
  # ────────────────────────────────────────────────────────────────────────────
  # DATA FILTERING AND VALIDATION
  # ────────────────────────────────────────────────────────────────────────────
  
  # Filter data for specific station, water mass, and variable combination
  data_subset <- data %>%
    filter(Station_Cluster == station, classification == water_mass, variable == !!variable)
  
  # Quality control: ensure minimum sample size for trend analysis
  if(nrow(data_subset) < 3) {
    cat("❌ Insufficient data points (n =", nrow(data_subset), 
        ") for Station", station, "| Water Mass:", water_mass, "| Variable:", variable, "\n")
    return(NULL)
  }
  
  # ────────────────────────────────────────────────────────────────────────────
  # OUTLIER DETECTION AND REMOVAL
  # ────────────────────────────────────────────────────────────────────────────
  
  # Apply 3×IQR rule to remove extreme outliers that may skew trend analysis
  Q1 <- quantile(data_subset$Value, 0.25, na.rm = TRUE)
  Q3 <- quantile(data_subset$Value, 0.75, na.rm = TRUE)
  data_subset <- filter(data_subset, 
                        Value >= Q1 - 3*(Q3-Q1), 
                        Value <= Q3 + 3*(Q3-Q1))
  
  # ────────────────────────────────────────────────────────────────────────────
  # STATISTICAL ANALYSIS
  # ────────────────────────────────────────────────────────────────────────────
  
  # Extract non-parametric trend statistics (Sen's slope and Mann-Kendall test)
  slope_sen <- first(data_subset$sen_slope)
  p_sen <- first(data_subset$p_value_mk)
  p_signif <- first(data_subset$p_value_signif)
  
  # Calculate Sen's slope intercept for trend line visualization
  intercept_sen <- median(data_subset$Value) - slope_sen * median(data_subset$dec_year)
  
  # Parametric comparison: Ordinary Least Squares regression
  mod <- lm(Value ~ dec_year, data = data_subset)
  slope_ols <- coef(mod)[2]
  p_ols <- summary(mod)$coeff[2,4]
  r2_ols <- summary(mod)$r.squared
  
  # ────────────────────────────────────────────────────────────────────────────
  # PLOT AESTHETICS AND LAYOUT PREPARATION
  # ────────────────────────────────────────────────────────────────────────────
  
  # Dynamic Y-axis scaling with 15% margin for visual clarity
  y_range <- range(data_subset$Value)
  y_margin <- diff(y_range) * 0.15
  var_name <- first(data_subset$variable)
  
  # Adaptive X-axis tick spacing based on temporal span
  year_range <- range(data_subset$dec_year, na.rm = TRUE)
  year_min <- floor(year_range[1])
  year_max <- ceiling(year_range[2])
  year_span <- year_max - year_min
  
  # Optimize tick frequency to prevent overcrowding
  if(year_span <= 10) {
    year_breaks <- seq(year_min, year_max, by = 1)     # Annual ticks for short series
  } else if(year_span <= 20) {
    year_breaks <- seq(year_min, year_max, by = 2)     # Biennial ticks for medium series
  } else {
    year_breaks <- seq(year_min, year_max, by = 5)     # 5-year ticks for long series
  }
  
  # ────────────────────────────────────────────────────────────────────────────
  # PUBLICATION-QUALITY PLOT GENERATION
  # ────────────────────────────────────────────────────────────────────────────
  
  ggplot(data_subset, aes(dec_year, Value)) +
    
    # Data points
    geom_point(size = 5) +
    
    # Non-parametric trend line (Sen's slope) - Primary trend indicator
    geom_abline(intercept = intercept_sen, slope = slope_sen, 
                color = "firebrick", linewidth = 1.2) +
    
    # Parametric trend line (OLS) - Comparison reference
    geom_smooth(method = "lm", se = F, color = "steelblue", 
                linetype = "dashed", alpha = 0.6) +
    
    # Statistical annotations - Sen's method (robust, non-parametric)
    annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.25,
             label = paste0("Sen: ", round(slope_sen, 4), ", p=", signif(p_sen, 3), " ", p_signif,", n=", nrow(data_subset)),
             color = "firebrick", size = 5, fontface = "bold") +
    
    # Statistical annotations - OLS method (parametric comparison)
    annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 3,
             label = paste0("OLS: ", round(slope_ols, 4), ", p=", signif(p_ols, 3), ", R²=", round(r2_ols, 3)),
             color = "steelblue", size = 5) +
    
    # Plot labels and titles
    labs(
      title = paste0("Station ", station, " | ", water_mass, " | ", var_labels[var_name]),
      x = "Year", 
      y = var_labels[var_name]
    ) +
    
    # X-axis configuration with adaptive ticks
    scale_x_continuous(
      breaks = year_breaks,
      labels = year_breaks,
      expand = expansion(mult = c(0.05, 0.05))  # Minimal margins
    ) +
    
    # Y-axis limits with calculated margins
    coord_cartesian(ylim = y_range + c(-y_margin, y_margin)) +
    
    # Professional theme for scientific publication
    theme(
      axis.text = element_text(size = 15, colour = "black"),
      axis.title = element_text(size = 15, colour = "black"),
      axis.line = element_line(linewidth = .75, colour = "black"),
      axis.ticks = element_line(linewidth = .75, colour = "black"),
      title = element_text(size = 20, colour = "black"),
      panel.background = element_blank(),
      panel.spacing = unit(2, "lines"),
      strip.text = element_text(size = 17.5, face = "bold")
    )
}

# ══════════════════════════════════════════════════════════════════════════════
# USAGE EXAMPLE
# ══════════════════════════════════════════════════════════════════════════════
# 
# # Generate trend plot for specific combination
# plot_trend(data = my_oceanographic_data, 
#            station = 113, 
#            water_mass = "PMW", 
#            variable = "Nitrate")
# 
# # Expected var_labels object (define globally):
# var_labels <- c(
#   "Temperature" = "Temperature (°C)",
#   "Salinity" = "Salinity (PSU)",
#   "Nitrate" = "Nitrate (µM)",
#   "Phosphate" = "Phosphate (µM)",
#   "Silicate" = "Silicate (µM)",
#   "Oxconc" = "Oxygen (µM)",
#   "Nstar" = "N*"
# )
# ══════════════════════════════════════════════════════════════════════════════
