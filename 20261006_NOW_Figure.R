# ══════════════════════════════════════════════════════════════════════════════
# WATER MASS VISUALIZATION AND CLASSIFICATION ANALYSIS
# ══════════════════════════════════════════════════════════════════════════════
#
# This section generates diagnostic plots for water mass identification and
# classification validation in the North Water (NOW) polynya system.
#
# Author: Nicolas Schiffrine
# Institution: Takuvik Joint International Laboratory, Université Laval
# ══════════════════════════════════════════════════════════════════════════════

# ──────────────────────────────────────────────────────────────────────────────
# DEPENDENCIES AND DATA LOADING
# ──────────────────────────────────────────────────────────────────────────────

# External datasets
NOW_Bathy <- read.csv("/Users/nicolas/Library/CloudStorage/OneDrive-UniversitéLaval/DataBase_Nut/Data/NOW_bathymetry.csv")

# Analysis functions and utilities
source("Script_R/sen_slope.R") # Mann-Kendall trend analysis with Sen's slope
source("Script_R/plot_section.R") # Oceanographic section plotting utilities
source("DataBase_Nut/20260209/20260209_PlotWaterMassSection.R") # σ₀-spicity water mass classification
source("DataBase_Nut/20260209/20260209_SensPlot.R") # Sensitivity analysis for Pacific Water tracers
source("DataBase_Nut/20260209/20260209_pw_regime_calc.R") # Pacific Water regime classification functions

# ──────────────────────────────────────────────────────────────────────────────
# OUTPUT CONFIGURATION
# ──────────────────────────────────────────────────────────────────────────────

file.format <- ".png" # Publication-ready format
file.path <- "DataBase_Nut/20260209/Figure/" # Output directory for all figures

# ──────────────────────────────────────────────────────────────────────────────
# BATHYMETRY DATA GENERATION (ARCHIVED CODE)
# ──────────────────────────────────────────────────────────────────────────────
#
# The following code generates bathymetry ribbon data for section plots.
# Currently archived as pre-computed bathymetry is loaded from CSV above.
# Uncomment and run if bathymetry data needs to be regenerated.
#
# Bathy <- tibble(
#   decimalLongitude = seq(
#     round(min(unique(NOW$decimalLongitude)), 1),
#     round(max(unique(NOW$decimalLongitude)), 1),
#     by = 0.05                                                             # High resolution for smooth bathymetry
#   ),
#   decimalLatitude = rep(mean(unique(NOW$decimalLatitude)))                # Constant latitude across transect
# ) %>%
#   marmap::get_depth() %>%                                                 # Query ETOPO1 global bathymetry
#   rename(Depth = depth)                                                   # Consistent naming convention
#
# write.csv(Bathy, "DataBase_Nut/Data/NOW_bathymetry.csv", row.names = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# GLOBAL PLOTTING CONFIGURATION
# ══════════════════════════════════════════════════════════════════════════════

# ──────────────────────────────────────────────────────────────────────────────
# WATER MASS COLOR PALETTE
# ──────────────────────────────────────────────────────────────────────────────

#' Colorblind-friendly water mass palette following Okabe-Ito scheme
#'
#' Ensures accessibility while maintaining oceanographic conventions.
#' Colors selected for maximum contrast and scientific clarity in publications.
#'
#' @details Water mass abbreviations and oceanographic context:
#' \itemize{
#'   \item PMW:  Polar Mode water Water (shallow, Pacific-origin)
#'   \item UHW:  Upper Halocline Water (intermediate depth)
#'   \item BBPW: Baffin Bay Polar Water (local formation)
#'   \item SPMW: Sub-Polar Mode Water (deeper, Atlantic-influenced)
#'   \item CBAW: Canada Basin Atlantic Water (deepest, Atlantic-origin)
#' }
#'
#' @reference Generated using: paletteer_d("basetheme::clean", n = 10)
wm_palette <- c(
  "PMW"  = "#5DA5DAFF", # Light blue - Pacific waters
  "UHW"  = "#FAA43AFF", # Orange - Halocline waters
  # "CBAW" = "#F15854FF", # Red - Deep Atlantic waters
  "BBPW" = "#60BD68FF", # Green - Baffin Bay waters
  "SPMW" = "#B276B2FF" # Purple - Sub-polar waters
)

# ──────────────────────────────────────────────────────────────────────────────
# OCEANOGRAPHIC ORDERING AND LABELING SYSTEMS
# ──────────────────────────────────────────────────────────────────────────────

#' Water mass depth stratification order (surface → deep)
#'
#' Reflects typical Arctic Ocean vertical structure in the NOW region.
#' Order based on characteristic depth ranges and density stratification.
Water_order <- c("PMW", "UHW", "BBPW", "SPMW") # , "CBAW")

#' Biogeochemical variable hierarchy for consistent plotting
#'
#' Groups physical properties first, then nutrients in order of
#' Biogeochemical variable hierarchy for consistent plotting
#'
#' Groups physical properties first, then nutrients in order of
#' oceanographic importance for Arctic biogeochemical cycling.
var_levels <- c(
  "CT", "SA", # Physical properties (TEOS-10)
  "Nitrate", "Phosphate", "Silicate", # "Oxconc", # Major nutrients + oxygen
  "NP", "SiN", "SiP", # Nutrient stoichiometric ratios
  "Nstar", "fPW" # Derived nutrient tracers
)

#' Publication-ready variable labels with proper units and formatting
#'
#' Follows TEOS-10 (IOC et al., 2010) conventions for axis labels and legends.
#' Uses expression() for proper subscript/superscript rendering.
var_labels <- list(
  # Physical properties (TEOS-10)
  CT = expression(Theta ~ "(°C)"), # Conservative Temperature
  SA = expression(S[A] ~ "(g" ~ kg^-1 * ")"), # Absolute Salinity

  # Inorganic nutrients and dissolved oxygen (gravimetric units)
  Nitrate = expression(NO[3]^"-" ~ "(" * mu * "mol" ~ kg^-1 * ")"),
  Nitrite = expression(NO[2]^"-" ~ "(" * mu * "mol" ~ kg^-1 * ")"),
  Ammonium = expression(NH[4]^"+" ~ "(" * mu * "mol" ~ kg^-1 * ")"),
  TIN = expression("TIN (" * mu * "mol" ~ kg^-1 * ")"),
  Phosphate = expression(PO[4]^"3-" ~ "(" * mu * "mol" ~ kg^-1 * ")"),
  Silicate = expression(Si(OH)[4] ~ "(" * mu * "mol" ~ kg^-1 * ")"),
  Oxconc = expression(O[2] ~ "(" * mu * "mol" ~ kg^-1 * ")"),

  # Nutrient ratios (molar, dimensionless)
  NP = "N:P",
  SiN = "Si:N",
  SiP = "Si:P",

  # Derived tracers
  Nstar = expression(N^"*" ~ "(" * mu * "mol" ~ kg^-1 * ")"), # Gruber & Sarmiento (1997)
  NO = expression("NO (" * mu * "mol" ~ kg^-1 * ")"), # 9·NO3 + O2 (Broecker, 1974)
  PO = expression("PO (" * mu * "mol" ~ kg^-1 * ")"), # 135·PO4 + O2
  fPW = expression(f[PW]) # Pacific Water fraction (0–1)
)

# ──────────────────────────────────────────────────────────────────────────────
# PLOTTING THEMES AND UTILITIES
# ──────────────────────────────────────────────────────────────────────────────

#' Base theme for multi-panel oceanographic figures
#'
#' Optimized for readability in scientific publications with multiple subplots.
#' Emphasizes data over decorative elements following Tufte principles.
base_theme <- theme(
  # Text elements sized for publication readability
  axis.text.x = element_text(size = 15, colour = "black", angle = 45, hjust = 1),
  axis.text.y = element_text(size = 15, colour = "black"),
  axis.title = element_text(size = 15, colour = "black"),
  strip.text = element_text(size = 17.5, face = "bold"),

  # Axis styling for professional appearance
  axis.line = element_line(linewidth = 0.75, colour = "black"),
  axis.ticks = element_line(linewidth = 0.75, colour = "black"),

  # Panel layout optimization
  panel.background = element_blank(), # Clean background
  panel.spacing = unit(2, "lines") # Adequate spacing between panels
)

#' Standard theme for individual oceanographic plots
#'
#' Based on theme_classic() with modifications for oceanographic conventions.
#' Removes unnecessary elements while maintaining professional appearance.
plot_theme <- theme_classic(base_size = 20) +
  theme(
    axis.text        = element_text(colour = "black"), # Black text for better contrast
    axis.ticks.x.top = element_blank(), # Remove top axis ticks
    axis.line.x.top  = element_blank() # Remove top axis line
  )

#' Empty plot placeholder for grid arrangements
#'
#' Used to fill unused positions in multi-panel layouts (e.g., 3×3 grids).
#' Maintains consistent spacing while avoiding visual clutter.
blank_plot <- ggplot() +
  geom_blank(aes(1, 1)) +
  theme_void()

# ──────────────────────────────────────────────────────────────────────────────
# STATISTICAL UTILITIES
# ──────────────────────────────────────────────────────────────────────────────

#' Convert p-values to significance symbols for scientific notation
#'
#' Follows standard conventions in oceanographic literature.
#' Uses conservative thresholds appropriate for environmental datasets.
#'
#' @param p Numeric vector of p-values
#' @return Character vector of significance symbols ("***", "**", "*", or NA)
#'
#' @details Significance levels:
#' \itemize{
#'   \item p < 0.001: *** (highly significant)
#'   \item p < 0.01:  **  (very significant)
#'   \item p < 0.05:  *   (significant)
#'   \item p ≥ 0.05:  NA  (not significant, blanked in plots)
#' }
# signif_p.value <- function(p) {
#   case_when(
#     p < 0.001 ~ "***",
#     p < 0.01  ~ "**",
#     p < 0.05  ~ "*",
#     TRUE      ~ NA_character_
#   )
# }

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 1: NORTH WATER TRANSECT OVERVIEW AND SAMPLING COVERAGE
# ══════════════════════════════════════════════════════════════════════════════
#
# This section generates a comprehensive overview figure showing the spatial
# and temporal distribution of oceanographic sampling across the NOW transect.
# Three complementary visualizations provide context for subsequent analyses:
# - Geographic map showing station locations
# - Temporal coverage showing sampling history (2005-2023)
# - Bathymetric context with depth distribution of measurements
#
# Output: Multi-panel figure suitable for manuscript Figure 1
# ══════════════════════════════════════════════════════════════════════════════

# ──────────────────────────────────────────────────────────────────────────────
# 1.1 GEOGRAPHIC MAP: NOW STATION LOCATIONS
# ──────────────────────────────────────────────────────────────────────────────

#' Create basemap of North Water polynya region
#'
#' Geographic context map showing station positions across the NOW transect.
#' Covers Baffin Bay to Nares Strait with optimal projection for Arctic regions.
#' Each station cluster is color-coded for identification across all panels.
#'
#' @details Geographic bounds optimized for NOW polynya system:
#' \itemize{
#'   \item Western boundary (-85°W): Captures full Baffin Bay extent
#'   \item Eastern boundary (-65°W): Includes Greenland coastal waters
#'   \item Southern boundary (75°N): Southern edge of polynya influence
#'   \item Northern boundary (77.5°N): Northern extent including Nares Strait
#' }
NOW_map <- basemap(
  limits = c(
    -85, # Western border: Full Baffin Bay coverage
    -65, # Eastern border: Greenland coastal margin
    75, # Southern border: Southern polynya extent
    77.5 # Northern border: Nares Strait region
  ),
  rotate = TRUE # Display longitudes as °W (oceanographic convention)
) +
  # Station positions with cluster-specific coloring
  ggspatial::geom_spatial_point(
    data = NOW,
    aes(
      x = decimalLongitude,
      y = decimalLatitude,
      fill = factor(Station_Cluster) # Cluster ID for consistent coloring across panels
    ),
    shape = 21, # Filled circles with black outline
    size = 5 # Large enough for publication visibility
  ) +
  # Colorblind-friendly station palette
  scale_fill_brewer(
    palette = "Set3", # Pastel colors for geographic clarity
    name = "Station"
  ) +
  # Professional cartographic styling
  theme(
    # Text sizing for publication readability
    axis.text = element_text(size = 17.5, colour = "black"),
    axis.title = element_text(size = 17.5, colour = "black"),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15),
    strip.text = element_text(size = 15, face = "bold"),

    # Professional map appearance
    panel.background = element_blank(), # Clean background for geographic clarity
    axis.ticks = element_line(linewidth = 1, colour = "black"),
    panel.border = element_rect(linewidth = 1.5, colour = "black"),
    legend.position = "right"
  )

# ──────────────────────────────────────────────────────────────────────────────
# 1.2 TEMPORAL COVERAGE: SAMPLING HISTORY ACROSS TRANSECT
# ──────────────────────────────────────────────────────────────────────────────

#' Temporal sampling distribution plot (Longitude × Time)
#'
#' Shows the evolution of sampling effort across the NOW transect from 2005-2023.
#' Reveals temporal biases, seasonal patterns, and data availability for trend analysis.
#' X-axis shows longitude gradient from Baffin Bay (west) to Nares Strait (east).
sampling_timeline <- ggplot(NOW, aes(x = -decimalLongitude, y = as.Date(eventDate))) +
  # Sampling points colored by station cluster
  geom_point(
    aes(fill = factor(Station_Cluster)),
    shape = 21, # Consistent with map symbols
    size = 5 # Publication-appropriate size
  ) +
  # Consistent color scheme across all panels
  scale_fill_brewer(
    palette = "Set3",
    name = "Station",
    na.value = "grey50" # Handle any missing cluster assignments
  ) +
  # Time axis formatting for multi-year dataset
  scale_y_date(
    date_breaks = "2 years", # Reasonable spacing for 18-year dataset
    date_labels = "%Y" # Year-only labels for clarity
  ) +
  # Longitude axis: reversed to match oceanographic convention (west → east)
  scale_x_reverse(
    breaks = seq(
      round(max(-NOW$decimalLongitude), 1) + 1.5, # Western extent with margin
      round(min(-NOW$decimalLongitude), 1) - 1.5, # Eastern extent with margin
      by = -1 # 1-degree intervals
    )
  ) +
  # Axis labels following oceanographic conventions
  labs(
    x = "Longitude (°W)",
    y = "Year"
  ) +
  # Consistent theme with map styling
  plot_theme +
  theme(
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 2),
    strip.background = element_rect(fill = "grey50"),
    strip.text       = element_text(face = "bold"),
    axis.line        = element_blank() # Clean borders instead of axis lines
  )

# ──────────────────────────────────────────────────────────────────────────────
# 1.3 BATHYMETRIC CONTEXT: DEPTH DISTRIBUTION WITH SEAFLOOR PROFILE
# ──────────────────────────────────────────────────────────────────────────────

#' Bathymetric section with sampling depth distribution
#'
#' Shows the relationship between sampling depths and seafloor bathymetry.
#' Critical for interpreting water mass distributions and benthic influences.
#' Grey ribbon represents seafloor, points show actual measurement depths.
bathymetric_section <- ggplot(NOW, aes(x = -decimalLongitude, y = Depth)) +
  # Sampling depths colored by station cluster
  geom_point(
    aes(fill = factor(Station_Cluster)),
    shape = 21,
    size = 2.5 # Slightly smaller to avoid overcrowding
  ) +
  # Seafloor bathymetry ribbon (grey background)
  geom_ribbon(
    data = NOW_Bathy,
    aes(
      x = -decimalLongitude,
      ymin = Depth, # Seafloor depth (negative values)
      ymax = rep(round(max(Depth), 1), length(Depth)) # Maximum depth for ribbon extent
    ),
    fill = "grey75", # Subtle background color
    inherit.aes = FALSE # Don't inherit point aesthetics
  ) +
  # Consistent station coloring
  scale_fill_brewer(palette = "Set3", name = "Station") +

  # Oceanographic depth convention: surface = 0, increasing downward
  scale_y_reverse() +

  # Longitude axis matching other panels
  scale_x_reverse(
    breaks = seq(
      round(max(-NOW$decimalLongitude), 1) + 1.5,
      round(min(-NOW$decimalLongitude), 1) - 1.5,
      by = -1
    )
  ) +
  # Standard oceanographic depth labeling
  labs(
    x = "Longitude (°W)",
    y = "Depth (m)"
  ) +
  # Consistent professional styling
  plot_theme +
  theme(
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 2),
    strip.background = element_rect(fill = "grey50"),
    strip.text       = element_text(face = "bold"),
    axis.line        = element_blank()
  )

# ──────────────────────────────────────────────────────────────────────────────
# 1.4 FIGURE COMPOSITION: MULTI-PANEL LAYOUT ASSEMBLY
# ──────────────────────────────────────────────────────────────────────────────

# get_legend <- function(myggplot) {
#   tmp <- ggplot_gtable(ggplot_build(myggplot))
#   leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
#   tmp$grobs[[leg]]
# }
#
# # Build a temporary plot to grab the legend (guide_legend with two columns)
# tmp <- Sampling +
#   guides(fill = guide_legend(ncol = 2))
#
# leg_tbl <- get_legend(tmp) # extract the legend as a grob
# leg_gg <- ggpubr::as_ggplot(leg_tbl) # turn it into a ggplot object

#' Remove individual legends for clean multi-panel layout
#'
#' Legends removed from individual panels to avoid redundancy.
#' Final figure uses unified legend or summary table in bottom-right panel.
panel_A <- NOW_map + theme(legend.position = "none")
panel_B <- sampling_timeline + theme(legend.position = "none")
panel_C <- bathymetric_section + theme(legend.position = "none")

#' Assemble publication-ready multi-panel overview figure
#'
#' Creates Figure 1 for manuscript: comprehensive NOW transect overview.
#' Layout optimized for two-column journal format with clear panel labeling.
#'
#' @details Panel arrangement:
#' \itemize{
#'   \item Panel A (top-left): Geographic map showing station locations
#'   \item Panel B (top-right): Temporal coverage (longitude × time)
#'   \item Panel C (bottom-left): Bathymetric section with sampling depths
#'   \item Panel D (bottom-right): Reserved for legend or summary table
#' }
NOW_overview <- ggpubr::ggarrange(
  panel_A, panel_B, panel_C, # Three main visualization panels
  # get_legend,
  labels = LETTERS, # Standard scientific panel labeling (A, B, C)
  nrow = 2, # Two-row layout for journal format
  ncol = 2, # Two-column layout
  widths = c(1, 1), # Equal width panels
  heights = c(1, 1), # Equal height panels
  align = "v" # Vertical alignment for professional appearance
)

# Display assembled figure for review
NOW_overview

#' Export high-resolution figure for manuscript submission
#'
#' Saves publication-ready figure with optimal dimensions and resolution.
#' Dimensions chosen for two-column journal format with clear readability.
# ggsave(
#   plot = NOW_overview,
#   filename = paste0("NOW_overview", file.format), # Uses global file format setting
#   path = file.path, # Uses global output directory
#   width = 15, # Width optimized for two-column format
#   height = 7.5, # Height maintaining good aspect ratio
#   dpi = 300 # Publication-quality resolution
# )

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 2: WATER MASS DIAGNOSTICS AND CLASSIFICATION VALIDATION
# ══════════════════════════════════════════════════════════════════════════════
#
# This section generates diagnostic plots for water mass identification and
# classification validation in the North Water (NOW) polynya system.
# Combines traditional T-S analysis with modern σ₀-spicity methodology
# following Huang et al. (2024) for robust Arctic water mass classification.
#
# Output: Multi-panel Figure 2 showing water mass diagnostic plots
# ══════════════════════════════════════════════════════════════════════════════

# ──────────────────────────────────────────────────────────────────────────────
# 2.1 TEMPERATURE-SALINITY DIAGRAM WITH ISOPYCNAL ANALYSIS
# ──────────────────────────────────────────────────────────────────────────────

#' Classical T-S diagram for water mass identification (Figure 2A)
#'
#' Traditional oceanographic approach using Temperature-Salinity relationships
#' with overlaid potential density (σ₀) contours for water mass characterization.
#' Depth-coded observations reveal vertical water column structure and mixing.
#'
#' @details Key diagnostic features:
#' \itemize{
#'   \item Isopycnal contours: Calculated using gsw package for accuracy
#'   \item Water mass endmembers: Literature values from Burger et al. (2023)
#'   \item Depth gradient: Ocean.deep color scale following conventions
#'   \item Reference coordinates: Median transect position for density calculations
#' }
#'
#' @note Arctic waters often exhibit complex T-S relationships due to:
#' - Sea ice formation/melting processes
#' - Strong haloclines from freshwater inputs
#' - Winter convection and summer stratification
TS_plot <- ggplot(NOW, aes(x = SA, y = CT, fill = Depth)) +
  # Main data points with depth gradient coloring
  geom_point(
    shape = 21, # Filled circles with outline
    size = 3.5, # Publication-appropriate size
    alpha = 0.75, # Slight transparency for overlapping points
    color = "grey25", # Dark outline for contrast
    stroke = 1 # Outline thickness
  ) +

  # Potential density isopycnal contours using accurate seawater calculations
  ggoce::geom_isopycnal(
    salinity_type    = "absolute", # Absolute Salinity (PSU)
    temperature_type = "conservative", # CT (not potential)
    ref_longitude    = median(NOW$decimalLongitude), # Reference longitude for calculations
    ref_latitude     = median(NOW$decimalLatitude), # Reference latitude for calculations
    trim_freezing    = FALSE # Include sub-freezing Arctic waters
  ) +

  # Water mass endmember positions from literature
  geom_point(
    data = wm_ref %>% filter(!watermass %in% "CBAW"),
    aes(x = SA, y = CT, shape = watermass),
    size = 5, # Prominent endmember symbols
    fill = "white", # White fill for visibility
    color = "black", # Black outline for contrast
    stroke = 1.2, # Thick outline for emphasis
    inherit.aes = FALSE # Don't inherit depth coloring
  ) +

  # Depth-based color scale using oceanographic conventions
  scale_fill_gradientn(
    colours = rev(ocean.deep(100)), # Deep ocean color palette (reversed)
    name = "Depth (m)" # Clear legend title
  ) +

  # Water mass endmember symbols
  scale_shape_manual(
    values = c(22, 23, 24, 25, 21), # Distinct geometric shapes
    name = "Water Mass Endmember"
  ) +

  # Axis labels following oceanographic conventions
  labs(
    x = expression(S[A] ~ "(g" ~ kg^-1 * ")"),
    y = expression(Theta ~ "(°C)")
  ) +

  # Legend customization for publication clarity
  guides(
    shape = guide_legend(
      override.aes = list(size = 4), # Smaller shapes in legend for space
      title = "Water Mass" # Concise legend title
    )
  ) +

  # Consistent theme across diagnostic plots
  plot_theme

# Display T-S diagnostic plot for review
TS_plot

# ──────────────────────────────────────────────────────────────────────────────
# 2.2 SIGMA-ZERO SPICITY CLASSIFICATION DIAGRAM
# ──────────────────────────────────────────────────────────────────────────────

#' σ₀-π₀ water mass classification plot (Figure 2B)
#'
#' Advanced water mass classification using σ₀-spicity methodology from
#' Huang et al. (2024). More robust than T-S analysis for Arctic systems
#' with complex mixing and strong stratification.
#'
#' @details Methodological advantages:
#' \itemize{
#'   \item σ₀-π₀ space: More conservative than T-S for mixing analysis
#'   \item Distance-based classification: Quantitative assignment to endmembers
#'   \item Arctic-specific: Optimized for polar water mass identification
#'   \item Mixing validation: Distance gradients show classification confidence
#' }
#'
#' @note Spicity (π₀) represents thermohaline variability orthogonal to density,
#' providing enhanced resolution of water mass mixing processes.

# Initialize base plot for multi-layer construction
spicity_base <- ggplot()

# Add each classified water mass with distance-based color gradients
for (mass_type in names(wm_palette)) {
  # Extract observations classified as current water mass
  mass_subset <- subset(NOW, classification == mass_type)

  # Add layer with independent color scale for each water mass
  spicity_base <- spicity_base +
    ggnewscale::new_scale_fill() + # Enable multiple fill scales
    geom_point(
      data = mass_subset,
      aes(x = sigma0, y = spicity, fill = distance), # Distance to nearest endmember
      size = 3.5, # Consistent point size
      shape = 21, # Filled circles
      color = "black", # Black outlines for contrast
      stroke = 1 # Outline thickness
    ) +
    # Distance-based color gradient from water mass color to white
    scale_fill_binned(
      low = wm_palette[mass_type], # Water mass-specific base color
      high = "white", # White for maximum distance
      limits = c(0, 0.6), # Standardized distance range
      name = paste(mass_type, "dist."), # Water mass-specific legend
      na.value = "grey50", # Grey for missing distances
      guide = guide_coloursteps(
        direction = "horizontal", # Horizontal color bars
        title.position = "top", # Title above color bar
        title.hjust = 0.5, # Center title
        barwidth = unit(3, "cm"), # Standard bar width
        barheight = unit(0.5, "cm") # Compact bar height
      )
    )
}

# Complete spicity plot with endmembers and mixed waters
spicity_plot <- spicity_base +
  # Add literature water mass endmember positions
  geom_point(
    data = wm_ref %>% filter(!watermass %in% "CBAW"),
    aes(x = sigma0, y = spicity, shape = watermass),
    size = 5, # Prominent endmember symbols
    stroke = 1.5, # Thick outlines
    color = "black", # Black outlines
    fill = "white", # White fill for visibility
    inherit.aes = FALSE # Independent aesthetics
  ) +

  # Add mixed/unclassified waters
  geom_point(
    data = subset(NOW, classification == "Mixed"),
    aes(x = sigma0, y = spicity),
    size = 3.5, # Standard point size
    shape = 21, # Filled circles
    fill = "grey", # Grey for mixed waters
    color = "black", # Black outlines
    stroke = 1, # Standard outline
    inherit.aes = FALSE # Independent aesthetics
  ) +

  # Endmember symbol mapping
  scale_shape_manual(
    values = c(22, 23, 24, 25, 21), # Distinct geometric shapes
    name = "Water Mass Endmember"
  ) +

  # Remove endmember legend (shown in Panel A)
  guides(shape = "none") +

  # Axis labels with proper scientific notation
  labs(
    x = expression("Potential Density" ~ sigma[0] ~ "(kg" ~ m^-3 * ")"), # σ₀ with proper subscript
    y = expression("Potential Density" ~ sigma[0] ~ "(kg" ~ m^-3 * ")") # π₀ with proper subscript
  ) +

  # Consistent theme with adjusted legend sizing
  plot_theme +
  theme(
    legend.text = element_text(size = 7.5), # Smaller text for multiple legends
    legend.title = element_text(size = 10) # Proportional title size
  )

# Display spicity diagnostic plot for review
print(spicity_plot)

# ──────────────────────────────────────────────────────────────────────────────
# 2.3 COMBINED WATER MASS DIAGNOSTIC FIGURE ASSEMBLY
# ──────────────────────────────────────────────────────────────────────────────

#' Assemble comprehensive water mass diagnostic figure (Figure 2)
#'
#' Combined T-S and σ₀-π₀ analysis providing complementary perspectives on
#' water mass structure in the NOW polynya system. Panel A shows classical
#' approach, Panel B shows modern classification methodology.
#'
#' @details Publication layout:
#' \itemize{
#'   \item Panel A: T-S diagram with isopycnals and endmembers
#'   \item Panel B: σ₀-π₀ classification with distance validation
#'   \item Horizontal layout: Optimal for two-column journal format
#'   \item Consistent scaling: Aligned axes for visual comparison
#' }
wm_plot <- ggpubr::ggarrange(
  TS_plot, # Panel A: Classical T-S analysis
  spicity_plot, # Panel B: Modern σ₀-π₀ classification
  nrow = 1, # Single row layout
  align = "hv", # Horizontal and vertical alignment
  labels = LETTERS # Standard panel labeling (A, B)
)

# Display combined diagnostic figure
# print(wm_plot)

# Export publication-ready combined figure
# ggsave(
#   plot = wm_plot,
#   filename = paste0("wm_plot", file.format), # Descriptive filename
#   path = file.path, # Global output directory
#   width = 22.5, # Wide format for side-by-side panels
#   height = 7.5, # Appropriate height for detail
#   dpi = 300 # Publication-quality resolution
# )

# ──────────────────────────────────────────────────────────────────────────────
# 2.4 CLASSIFICATION VALIDATION: DISTANCE DISTRIBUTION ANALYSIS
# ──────────────────────────────────────────────────────────────────────────────

#' Validate water mass classification through distance distribution analysis
#'
#' Examines the distribution of σ₀-π₀ distances from observations to nearest
#' water mass endmembers. Small distances indicate confident classification;
#' large distances suggest mixed waters or classification uncertainty.
#'
#' @details Validation metrics:
#' \itemize{
#'   \item Distance distribution: Histogram showing classification confidence
#'   \item Sample sizes: Number of observations per water mass
#'   \item Median distances: Central tendency of classification quality
#'   \item Faceted display: Individual assessment for each water mass
#' }

# Calculate comprehensive statistics for each classified water mass
wm_dist_stats <- NOW %>%
  filter(!classification %in% "Mixed") %>% # Exclude mixed waters from validation
  group_by(classification) %>% # Group by water mass type
  rstatix::get_summary_stats(distance, type = "full") %>% # Full statistical summary
  mutate(
    # Create annotation labels with key metrics
    stats_label = paste0("n = ", n, "; median = ", round(median, 3)) # Sample size and median distance
  )

# View(
#   wm_dist_stats %>%
#   mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
#   select(classification,n,median,iqr)
#   )

# EXPORT — WATER MASS CLASSIFICATION DISTANCE STATISTICS (Excel)
#' Summary table of classification distances per water mass
#'
#' Distance to the end-member reference used for water mass validation.
#' Values rounded to 3 decimals (consistent with the plot annotation).
wm_dist_summary <-
  NOW %>%
  group_by(classification) %>% # Group by water mass type
  rstatix::get_summary_stats(distance, show = c("n", "median", "iqr")) %>% # Keep reporting columns
  arrange(factor(classification, levels = Water_order)) # Canonical water mass order

View(wm_dist_summary) # Visual check

# ── Single-table export ──────────────────────────────────────────────────────
# writexl::write_xlsx(
#   x = wm_dist_summary, # Distance summary table
#   path = file.path(
#     file.path, # Global output directory
#     paste0("wm_distance_stats", ".xlsx")
#   ) # Descriptive filename
# )

# Create faceted validation plot
wm_dist_plot <- ggplot(
  data = NOW %>%
    filter(!classification %in% "Mixed"), # Focus on classified observations
  aes(x = distance, fill = classification)
) +
  # Distance distribution histograms
  geom_histogram(
    color = "white", # White bin outlines for clarity
    alpha = 0.8, # Slight transparency
    bins = 20 # Reasonable bin count for detail
  ) +

  # Separate panel for each water mass
  ggh4x::facet_wrap2(~classification, axes = "all", strip = strip_vanilla(clip = "on")) +

  # Water mass-specific coloring consistent with diagnostic plots
  scale_fill_manual(values = wm_palette) +

  # Add statistical annotations to each panel
  geom_text(
    data = wm_dist_stats,
    aes(x = Inf, y = Inf, label = stats_label),
    inherit.aes = FALSE, # Independent positioning
    hjust = 1.1, # Right-justified
    vjust = 1.5, # Top-positioned
    size = 5, # Readable annotation size
    fontface = "bold" # Emphasis for key statistics
  ) +

  # Clear axis labels for interpretation
  labs(
    x = "Distance to Closest Water Mass Endmember", # Descriptive x-axis
    y = "Count", # Clear y-axis meaning
    fill = "Water Mass" # Legend title
  ) +

  # Professional theme with facet customization
  plot_theme +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 2),
    axis.line = element_blank(),
    strip.background = element_rect(fill = "grey50", linewidth = 2), # Prominent facet headers
    strip.text = element_text(face = "bold"), # Bold water mass names
    legend.position = "none" # Remove redundant legend
  )

# Display validation plot for review
wm_dist_plot

# Export classification validation figure
# ggsave(
#   plot = wm_dist_plot,
#   filename = paste0("wm_dist", file.format), # Supplementary figure
#   path = file.path, # Global output directory
#   width = 15, # Wide format for multiple facets
#   height = 7.5, # Standard height
#   dpi = 300 # Publication-quality resolution
# )

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 3: VERTICAL STRUCTURE AND PROPERTY DISTRIBUTIONS
# ══════════════════════════════════════════════════════════════════════════════
#
# This section analyzes the vertical structure of physical and biogeochemical
# properties across the NOW transect. Examines how temperature, salinity,
# nutrients, and derived tracers vary by station and water mass classification.
#
# Key analyses:
# - Cross-transect gradients in water mass properties
# - Depth distribution of classified water masses
# - Statistical summary of property ranges within each water mass
#
# Output: Multi-panel Figure 3 showing vertical property structure
# ══════════════════════════════════════════════════════════════════════════════

# ──────────────────────────────────────────────────────────────────────────────
# 3.1 DATA PREPARATION AND STATISTICAL AGGREGATION
# ──────────────────────────────────────────────────────────────────────────────

#' Aggregate oceanographic properties by station and water mass classification
#'
#' Computes comprehensive statistical summaries for all measured variables
#' within each station × water mass combination. Excludes mixed waters to
#' focus on well-defined water mass properties.
#'
#' @details Statistical measures include:
#' \itemize{
#'   \item Sample size (n): Number of observations per group
#'   \item Central tendency: Median (robust to outliers)
#'   \item Variability: Q1, Q3 for interquartile range
#'   \item Extremes: Min, max for full data range
#' }
#'
#' @note Median and IQR chosen over mean and SD due to:
#' - Non-normal distributions common in oceanographic data
#' - Robustness to outliers from sensor drift or mixing events
#' - Better representation of water mass "core" properties
wm_stat <- NOW %>%
  filter(!classification %in% "Mixed") %>% # Exclude mixed/transitional waters
  group_by(Station_Cluster, classification) %>% # Group by station and water mass
  rstatix::get_summary_stats(any_of(var_levels), type = "full") %>% # Full statistical summary
  mutate(
    variable = factor(variable, levels = var_levels), # Maintain canonical variable order
    Var.plot = var_labels[as.character(variable)] # Apply readable plot labels
  ) %>%
  left_join(
    .,
    NOW %>% group_by(Station_Cluster) %>% summarise(lon_med = median(decimalLongitude, na.rm = TRUE), .groups = "drop")
  )

# EXPORT — WATER MASS SUMMARY STATISTICS (Excel)
#' Summary table of water-mass properties per station cluster
#'
#' Numeric values rounded to 2 decimals; central tendency and dispersion
#' formatted as "mean ± sd" and "median [q1 – q3]".
# ── Summary table (unchanged) ────────────────────────────────────────────────
wm_summary <- wm_stat %>%
  mutate(
    median_iq = paste0(median, " [", q1, " – ", q3, "]"),
    iq = paste0(" [", q1, " – ", q3, "]"),
    range_min_max     = paste0("[", min, " – ", max, "]")
  ) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
  select(
    Station_Cluster, classification, variable,
    n, median, iqr, median_iq,iq, range_min_max
  ) %>%
  arrange(variable, Station_Cluster, classification)

# ── One data frame per parameter ─────────────────────────────────────────────
wm_summary_list <- wm_summary %>%
  group_by(variable) %>%
  group_split() %>% # One tibble per variable
  set_names(map_chr(., ~ as.character(unique(.x$variable)))) %>% # Sheet name = variable
  map(~ select(.x, -variable)) # Drop redundant column

# ── Sheet 2: global statistics per water mass (all stations pooled) ──────────
wm_summary_global <- NOW %>%
  filter(!classification %in% "Mixed") %>% # Same exclusion as wm_stat
  group_by(classification) %>% # Pool all station clusters
  rstatix::get_summary_stats(any_of(var_levels), type = "full") %>%
  mutate(
    variable = factor(variable, levels = var_levels),
    Var.plot = var_labels[as.character(variable)]
  ) %>%
  ungroup() %>%
  mutate(
    median_iq = paste0(median, " [", q1, " – ", q3, "]"),
    iq = paste0(" [", q1, " – ", q3, "]"),
    range_min_max     = paste0("[", min, " – ", max, "]")
  ) %>%
  select(
    classification, variable,
    n, median, iqr, median_iq,iq, range_min_max
  ) %>%
  arrange(variable, classification)

# ── Optional: add the global-by-water-mass sheet first ───────────────────────
sheets <- c(list("Global_by_watermass" = wm_summary_global), wm_summary_list)

# ── Multi-sheet export ───────────────────────────────────────────────────────
# writexl::write_xlsx(
#   x    = sheets,
#   path = file.path(file.path, "wm_summary_stats.xlsx")
# )

# Table de correspondance station <-> longitude médiane
station_axis <- NOW %>%
  filter(!classification %in% "Mixed") %>% # Focus on classified observations
  group_by(Station_Cluster, classification) %>% # Group by station and water mass
  get_summary_stats(Depth, type = "full") %>% # Full depth statistics
  left_join(
    .,
    NOW %>% group_by(Station_Cluster) %>% summarise(lon_med = median(decimalLongitude, na.rm = TRUE), .groups = "drop")
  ) %>%
  distinct(Station_Cluster, lon_med) %>%
  arrange(lon_med)

# View(
#   wm_stat %>%
#   mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
#   select(Station_Cluster,classification,variable,n,median,iqr)
#   )

# ──────────────────────────────────────────────────────────────────────────────
# 3.2 VISUALIZATION FUNCTIONS FOR VERTICAL STRUCTURE ANALYSIS
# ──────────────────────────────────────────────────────────────────────────────

#' Create median ± IQR plots for oceanographic properties
#'
#' Generates publication-quality plots showing property distributions across
#' the NOW transect. Each water mass is color-coded and shape-coded for
#' clear identification. Uses robust statistics (median, IQR) appropriate
#' for oceanographic datasets with potential outliers.
#'
#' @param data Filtered dataset for single variable
#' @param variable Character string of variable name (for internal consistency)
#' @return ggplot2 object with station-by-station property comparison
#'
#' @details Design principles:
#' \itemize{
#'   \item Line ranges: Show IQR (25th-75th percentile) as measure of variability
#'   \item Point symbols: Show median values with water mass-specific shapes
#'   \item Position dodge: Separate overlapping water masses at same station
#'   \item Consistent coloring: Matches water mass diagnostic plots
#' }
plot_linerange <- function(data, show_lon_axis = FALSE) {
  station_axis <- data |>
    dplyr::distinct(Station_Cluster, lon_med) |>
    dplyr::arrange(lon_med) # west (most negative) → east

  p <- ggplot(data) +
    geom_linerange(
      aes(lon_med, median, # ← retiré le signe négatif
        ymin = q1, ymax = q3,
        colour = classification
      ),
      position = position_dodge2(width = 0.25, preserve = "single"),
      linewidth = 1.2
    ) +
    geom_point(
      aes(lon_med, median, # ← idem
        shape = classification,
        fill  = classification
      ),
      position = position_dodge2(width = 0.25, preserve = "single"),
      size = 3.5
    ) +
    scale_fill_manual(values = wm_palette, name = "Water Mass") +
    scale_color_manual(values = wm_palette, name = "Water Mass") +
    scale_shape_manual(name = "Water Mass", values = c(21, 22, 23, 24, 25)) +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
    labs(
      x = NULL,
      y = var_labels[[as.character(data$variable[1])]]
    ) +
    theme_classic(base_size = 15) +
    theme(
      axis.text = element_text(colour = "black"), # Black text for better contrast
      axis.ticks.x.top = element_blank(), # Remove top axis ticks
      axis.line.x.top = element_blank(), # Remove top axis line
      panel.border = element_rect(color = "black", fill = NA, linewidth = 2),
      axis.line = element_blank(),
      strip.background = element_rect(fill = "grey50", linewidth = 2), # Prominent facet headers
      strip.text = element_text(face = "bold"), # Bold water mass names
      legend.position = "none" # Remove redundant legend
    )

  if (show_lon_axis) {
    # Breaks cohérents avec les autres figures (section verticale)
    lon_breaks <- seq(
      floor(min(NOW$decimalLongitude)) - 1.5, # ex. -79
      ceiling(max(NOW$decimalLongitude)) + 1.5, # ex. -71
      by = 1
    )

    p <- p +
      scale_x_continuous(
        breaks = lon_breaks,
        labels = function(x) abs(round(x, 0)) # affiche |lon| → 78, 76, ...
      ) +
      guides(x = compose_stack(
        guide_axis_base(
          key = key_manual(
            aesthetic = station_axis$lon_med, # ← sans le -
            label     = station_axis$Station_Cluster
          ),
          theme = theme(
            axis.line.x = element_blank(),
            axis.text.x = element_text(size = 15, colour = "black", angle = 45, hjust = 1)
          )
        ),
        guide_axis_base(
          theme = theme(
            axis.line.x = element_line(linewidth = 0.75),
            axis.text.x = element_text(size = 15, colour = "black")
          )
        ),
        side.titles = c("Station", "Longitude (°W)")
      )) +
      theme_classic(base_size = 15) +
      theme(
        axis.text = element_text(colour = "black"), # Black text for better contrast
        axis.ticks.x.top = element_blank(), # Remove top axis ticks
        axis.line.x.top = element_blank(), # Remove top axis line
        panel.border = element_rect(color = "black", fill = NA, linewidth = 2),
        axis.line = element_blank(),
        strip.background = element_rect(fill = "grey50", linewidth = 2), # Prominent facet headers
        strip.text = element_text(face = "bold"), # Bold water mass names
        legend.position = "none", # Remove redundant legend
        plot.margin = margin(t = 5, r = 10, b = 5, l = 35, unit = "pt")
      )
  } else {
    p <- p +
      scale_x_continuous(
        breaks = station_axis$lon_med, # ← sans le -
        labels = station_axis$Station_Cluster
      ) +
      labs(x = NULL) +
      theme_classic(base_size = 15) +
      theme(
        axis.text.x = element_text(size = 15, colour = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 15, colour = "black"), # Black text for better contrast
        axis.ticks.x.top = element_blank(), # Remove top axis ticks
        axis.line.x.top = element_blank(), # Remove top axis line
        panel.border = element_rect(color = "black", fill = NA, linewidth = 2),
        axis.line = element_blank(),
        strip.background = element_rect(fill = "grey50", linewidth = 2), # Prominent facet headers
        strip.text = element_text(face = "bold"), # Bold water mass names
        legend.position = "none" # Remove redundant legend
      )
  }

  p
}

#' Assemble multi-panel figure grid with shared legend for water mass analysis
#'
#' Creates comprehensive figure showing all oceanographic properties
#' in organized grid layout. Optimized for manuscript Figure 3 with
#' professional spacing and labeling to reveal cross-transect variability
#' within classified water masses.
#'
#' @param plot_function Function to generate individual property panels
#' @param dataset Statistical summary dataset containing water mass properties
#' @param layout_config List containing variables, positions, and grid dimensions
#' @return Combined ggplot2 object with multiple panels and shared legend
#'
#' @details Grid organization optimized for oceanographic manuscript:
#' \itemize{
#'   \item Row 1: Physical properties (CT, SA) + strategic blanks
#'   \item Row 2: Primary nutrients (Nitrate, Phosphate, Silicate) revealing biogeochemical gradients
#'   \item Row 3: Derived tracers (N*) indicating water mass provenance + balanced layout
#'   \item Right side: Shared legend showing water mass classifications
#'   \item Panel labels: Sequential lettering (A-F) for manuscript referencing
#' }
#'
#' @examples
#' # Define layout for 6 variables in 3×3 grid
#' layout_config <- list(
#'   variables = c("CT", "SA", "Nitrate", "Phosphate", "Silicate", "Nstar"),
#'   positions = c(1, 2, 4, 5, 6, 7),
#'   nrow = 3, ncol = 3
#' )
#' grid <- make_grid(plot_linerange, wm_stat, layout_config)
#'
make_grid <- function(plot_function, dataset, layout_config) {
  blank_plot <- ggplot() +
    theme_void()

  variables <- layout_config$variables
  positions <- layout_config$positions
  nrow <- layout_config$nrow
  ncol <- layout_config$ncol

  total_slots <- nrow * ncol
  plots_list <- rep(list(blank_plot), total_slots)
  labels <- rep("", total_slots)

  # Identify bottom panels per column: for each column, find the LAST
  # occupied position (highest row index among filled slots)
  col_of <- function(p) ((p - 1) %% ncol) + 1
  row_of <- function(p) ((p - 1) %/% ncol) + 1

  occupied <- positions
  bottom_pos <- sapply(1:ncol, function(cc) {
    in_col <- occupied[col_of(occupied) == cc]
    if (length(in_col) == 0) NA_integer_ else in_col[which.max(row_of(in_col))]
  })
  bottom_pos <- bottom_pos[!is.na(bottom_pos)]

  for (i in seq_along(variables)) {
    pos <- positions[i]
    if (!is.na(pos) && pos <= total_slots) {
      plots_list[[pos]] <- plot_function(
        filter(dataset, variable == variables[i]),
        show_lon_axis = pos %in% bottom_pos
      )
      labels[pos] <- LETTERS[i]
    }
  }

  do.call(ggarrange, c(
    plots_list,
    list(
      labels        = labels,
      nrow          = nrow,
      ncol          = ncol,
      align         = "hv",
      common.legend = TRUE,
      legend        = "right"
    )
  ))
}

# ──────────────────────────────────────────────────────────────────────────────
# 3.3 GENERATE VERTICAL STRUCTURE ANALYSIS FIGURE (MANUSCRIPT FIGURE 3)
# ──────────────────────────────────────────────────────────────────────────────

#' Create comprehensive vertical structure figure for NOW transect analysis
#'
#' Multi-panel visualization revealing how physical and biogeochemical properties
#' vary across the North Water (NOW) transect within each classified water mass.
#'
#' Scientific objectives:
#' - Document cross-transect gradients in water mass properties
#' - Validate consistency of water mass classification across stations
#' - Identify biogeochemical signatures distinguishing water masses
#' - Reveal spatial patterns in nutrient distributions and water mass provenance
#'
#' Figure design rationale:
#' - Physical properties (T, S) establish water mass framework
#' - Nutrients (NO₃, PO₄, SiO₄) reveal biogeochemical gradients
#' - N* tracer indicates Pacific vs Atlantic water mass influence
#' - 3×3 grid provides balanced layout with strategic spacing

# Define layout configuration for optimal scientific presentation
layout_config <- list(
  # Variables ordered by oceanographic importance and logical grouping
  variables = c(
    "SA", "CT",
    "Nitrate", "Phosphate", "Silicate",
    # "NP","SiN","SiP",
    "Nstar", "fPW"
  ),

  # Positions create thematic grouping:
  # Row 1: Physical properties (1,2) + blank (3)
  # Row 2: Primary nutrients (4,5,6)
  # Row 3: Primary nutrients (7,8,9)
  # Row 4: Derived tracer (10,11) + blanks (12)
  positions = c(
    1, 2,
    4, 5, 6,
    7, 8
  ),

  # Grid dimensions optimized for manuscript format
  nrow = 3, # Three rows for thematic organization
  ncol = 3 # Three columns for balanced appearance
)

# Generate comprehensive water mass vertical structure figure
# Uses plot_linerange to show median ± IQR for each water mass at each station
wm_vertical_structure <- make_grid(
  plot_function = plot_linerange, # Function generating individual property panels
  dataset = wm_stat, # Statistical summary of water mass properties
  layout_config = layout_config # Layout configuration defined above
)

# Display the assembled multi-panel figure
# This reveals cross-transect patterns in water mass properties and
# validates the consistency of our water mass classification scheme
wm_vertical_structure

# Export publication-ready figure with optimal specifications
# Uncomment when ready for manuscript submission
# ggsave(
#   plot = wm_vertical_structure,
#   filename = paste0("wm_vertical_structure", file.format), # Descriptive filename for manuscript
#   path = file.path, # Global output directory
#   width = 20, # Wide format accommodates multi-panel layout and legend
#   height = 10, # Sufficient height for three-row scientific presentation
#   dpi = 300 # Publication-quality resolution for journal submission
# )

# ──────────────────────────────────────────────────────────────────────────────
# 3.4 WATER MASS DEPTH DISTRIBUTION ANALYSIS
# ──────────────────────────────────────────────────────────────────────────────

#' Analyze depth distribution of classified water masses
#'
#' Examines the vertical positioning of water masses across the NOW transect.
#' Critical for validating oceanographic understanding of water column
#' stratification and identifying potential misclassifications.
#'
#' @details Expected depth patterns:
#' \itemize{
#'   \item PMW (Pacific Mode Water): Surface to ~100m
#'   \item UHW (Upper Halocline Water): ~100-300m
#'   \item BBPW (Baffin Bay Polar Water): ~200-500m
#'   \item SPMW (Subpolar Mode Water): >400m
#'   \item CBAW (Canada Basin Atlantic Water): >500m (deepest)
#' }

# ── Station-level statistics (with median longitude for ordering) ────────────
wm_depth_stats <- NOW %>%
  filter(!classification %in% "Mixed") %>%
  group_by(Station_Cluster, classification) %>%
  get_summary_stats(Depth, type = "full") %>%
  ungroup() %>%
  left_join(
    NOW %>%
      group_by(Station_Cluster) %>%
      summarise(lon_med = median(decimalLongitude, na.rm = TRUE), .groups = "drop"),
    by = "Station_Cluster"
  )

# ── Formatting helper (no mean / SD) ─────────────────────────────────────────
fmt_depth <- function(df) {
  df %>%
    mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
    mutate(
      median_iq = paste0(median, " [", q1, " – ", q3, "]"),   # Median [IQR]
      range     = paste0("[", min, " – ", max, "]")           # Observed range
    )
}

# ── Sheet 1: global statistics per water mass (all stations pooled) ──────────
wm_depth_global <- NOW %>%
  filter(!classification %in% "Mixed") %>%
  group_by(classification) %>%
  get_summary_stats(Depth, type = "full") %>%
  ungroup() %>%
  fmt_depth() %>%
  select(classification, n, median, iqr, median_iq, range) %>%
  arrange(classification)

# ── Sheet 2: detail by station cluster × water mass ──────────────────────────
wm_depth_detail <- wm_depth_stats %>%
  fmt_depth() %>%
  arrange(lon_med, Station_Cluster, classification) %>%
  select(Station_Cluster, classification, n, median_iq, range)

# ── Export ───────────────────────────────────────────────────────────────────
# writexl::write_xlsx(
#   x = list(
#     "Global_by_watermass" = wm_depth_global,
#     "By_station"          = wm_depth_detail
#   ),
#   path = file.path(file.path, "wm_depth_stats.xlsx")
# )

# Table de correspondance station <-> longitude médiane
station_axis <- wm_depth_stats %>%
  distinct(Station_Cluster, lon_med) %>%
  arrange(lon_med)

# Create water mass depth distribution plot
wm_depth_plot <- ggplot(wm_depth_stats) +
  # Depth range as vertical line segments (IQR)
  geom_linerange(
    aes(
      -lon_med, median, # X: longitude, Y: median depth
      ymin = q1, ymax = q3, # Range: depth IQR
      colour = classification # Color by water mass
    ),
    # position  = position_dodge2(width = 0.25, preserve = "single"),
    linewidth = 1.2
  ) +

  # Median depths as filled points
  geom_point(
    aes(
      -lon_med, median,
      fill = classification,
      shape = classification
    ),
    # position = position_dodge2(width = 0.25, preserve = "single"),
    size = 3.5
  ) +

  # Consistent water mass aesthetics
  scale_fill_manual(values = wm_palette, name = "Water Mass") +
  scale_color_manual(values = wm_palette, name = "Water Mass") +
  scale_shape_manual(
    name = "Water Mass",
    values = c(21, 22, 23, 24, 25) # Water mass-specific shapes
  ) +
  scale_x_reverse(
    breaks = seq(
      round(max(-NOW$decimalLongitude), 1) + 1.5,
      round(min(-NOW$decimalLongitude), 1) - 1.5,
      by = -1
    )
  ) +

  # Inverted y-axis following oceanographic convention (depth increases downward)
  scale_y_reverse(expand = expansion(mult = c(0.25, 0.25))) +

  # Double axe X empilé : Longitude puis Station, titres à droite
  guides(x = compose_stack(
    guide_axis_base(
      key = key_manual(
        aesthetic = -station_axis$lon_med,
        label     = station_axis$Station_Cluster
      ),
      theme = theme(axis.line.x = element_blank())
    ),
    guide_axis_base(
      theme = theme(axis.line.x = element_line(linewidth = 0.75))
    ),
    side.titles = c("Station", "Longitude (°W)")
  )) +

  # Axis labels (y seulement ; x géré par les guides)
  labs(
    x = NULL,
    y = "Depth (m)"
  ) +

  # Enhanced theme for standalone figure
  plot_theme +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 2), # Prominent border
    strip.background = element_rect(fill = "grey50"), # Grey strip background
    strip.text = element_text(face = "bold"), # Bold strip text
    plot.margin = margin_auto(3, unit = "cm"),
    axis.line.x = element_line(color = "black", linewidth = 0.6),
    axis.line.y = element_blank()
  )

# Display water mass depth structure
wm_depth_plot

# Export depth structure figure
# ggsave(
#   plot = wm_depth_plot,
#   filename = paste0("wm_depth", file.format), # Supplementary figure
#   path = file.path, # Global output directory
#   width = 15, # Standard width for single-panel
#   height = 7.5, # Appropriate height for depth range
#   dpi = 300 # Publication-quality resolution
# )

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 4: TEMPORAL TREND ANALYSIS AND DECADAL VARIABILITY
# ══════════════════════════════════════════════════════════════════════════════
#
# This section quantifies multi-decadal trends (2005-2023) in physical and
# biogeochemical properties across the NOW transect using robust non-parametric
# statistical methods. Employs Mann-Kendall trend detection with Sen's slope
# estimation for trend magnitude quantification.
#
# Key analyses:
# - Station-by-station trend detection for all water masses
# - Magnitude estimation using Sen's slope (robust to outliers)
# - Significance testing with Mann-Kendall (no distributional assumptions)
# - Spatial patterns in temporal change across the transect
#
# Output: Multi-panel Figure 4 showing temporal trend heatmaps
# ══════════════════════════════════════════════════════════════════════════════

# ──────────────────────────────────────────────────────────────────────────────
# 4.1 TEMPORAL TREND CALCULATION AND STATISTICAL TESTING
# ──────────────────────────────────────────────────────────────────────────────

#' Calculate decadal trends using robust non-parametric methods
#'
#' Applies Mann-Kendall trend test combined with Sen's slope estimation to
#' detect and quantify temporal changes in oceanographic properties. This
#' approach is particularly suitable for environmental time series with:
#' - Non-normal distributions common in oceanographic data
#' - Potential outliers from extreme events or sensor issues
#' - Irregular sampling intervals across years and stations
#' - Small sample sizes in some station×water mass combinations
#'
#' @details Statistical methodology:
#' \itemize{
#'   \item Mann-Kendall test: Non-parametric trend detection (no normality assumption)
#'   \item Sen's slope: Robust trend magnitude (median of pairwise slopes)
#'   \item Minimum sample size: n≥3 for reliable trend estimation
#'   \item Significance threshold: p<0.05 with multiple testing awareness
#' }
#'
#' @note Arctic oceanography considerations:
#' - Sea ice variability affects sampling frequency and bias
#' - Climate oscillations (NAO, AO) may create non-linear trends
#' - Water mass source changes require station-specific analysis
temporal_trends_raw <- NOW %>%
  # Prepare data for trend analysis
  select(Station_Cluster, eventDate, classification, all_of(var_levels)) %>%
  pivot_longer(any_of(var_levels), names_to = "variable", values_to = "value") %>%
  drop_na(value) %>% # Remove missing observations

  # Group by analysis units
  group_by(Station_Cluster, classification, variable) %>%
  filter(n() >= 3) %>% # Minimum sample size for trends

  # Apply Mann-Kendall test and Sen's slope to each group
  group_modify(~ sen_slope(.x, eventDate, value)) %>% # Custom function for trend analysis
  ungroup()

# ──────────────────────────────────────────────────────────────────────────────
# 4.2 TREND DATA PROCESSING AND SIGNIFICANCE ASSESSMENT
# ──────────────────────────────────────────────────────────────────────────────

#' Process trend results for visualization and interpretation
#'
#' Cleans trend calculation outputs and applies significance filtering for
#' publication-quality visualization. Non-significant trends are masked to
#' emphasize meaningful temporal changes while avoiding false discoveries.
#'
#' @details Processing steps:
#' \itemize{
#'   \item Infinite values: Replace with NA (numerical stability)
#'   \item Significance symbols: Standard notation (***, **, *, ns)
#'   \item Masking: Hide non-significant trends in heatmap visualization
#'   \item Standardization: Consistent variable ordering and labeling
#' }
temporal_trends_processed <- temporal_trends_raw %>%
  # Clean numerical artifacts from trend calculations
  mutate(across(where(is.numeric), ~ if_else(is.infinite(.x), NA_real_, .x))) %>%
  # Generate significance symbols and plotting variables
  mutate(
    # Standard significance notation for scientific publication
    p_value_signif = case_when(
      is.na(sen_slope) ~ NA_character_,
      p_value_mk < 0.001 ~ "***", # Highly significant
      p_value_mk < 0.01 ~ "**", # Very significant
      p_value_mk < 0.05 ~ "*", # Significant
      TRUE ~ NA_character_ # Non-significant (masked)
    ),

    # Plotting variables (NA for non-significant trends)
    sen_slope_plot = if_else(is.na(p_value_signif), NA_real_, sen_slope), # Slope for significant trends only
    z_mk_plot = if_else(is.na(p_value_signif), NA_real_, z_mk), # Test statistic for significant trends

    # Standardize factor levels for consistent plotting
    variable = factor(variable, levels = var_levels), # Canonical variable order
    Var.plot = factor(var_labels[as.character(variable)], levels = var_labels) # Readable labels
  )


# EXPORT — TEMPORAL TREND SUMMARY TABLE
#' Summary table of temporal trends (Sen's slope, 95% CI, significance)
#'
#' Numeric values rounded to 2 decimals; CI formatted as "[lower – upper]".
trend_summary <-
  temporal_trends_processed %>%
  mutate(across(where(is.numeric), ~ round(.x, 2))) %>% # Round numeric columns
  mutate(ci = paste0("[", ci_lower, " – ", ci_upper, "]")) %>% # Format confidence interval
  select(
    Station_Cluster, classification, variable,
    n_obs, sen_slope, ci, p_value_signif
  ) # Keep reporting columns

View(trend_summary) # Visual check

# writexl::write_xlsx(
#   x = trend_summary, # Trend summary table
#   path = file.path(
#     file.path, # Global output directory
#     paste0("temporal_trends_summary", ".xlsx")
#   ) # Descriptive filename
# )

# ──────────────────────────────────────────────────────────────────────────────
# 4.3 TREND VISUALIZATION FUNCTIONS
# ──────────────────────────────────────────────────────────────────────────────

#' Create heatmap tiles for trend visualization
#'
#' Generates publication-quality heatmap showing trend magnitudes across
#' stations and water masses. Color intensity represents Sen's slope magnitude,
#' while overlaid symbols indicate statistical significance levels.
#'
#' @param data Filtered dataset for single variable
#' @param variable Variable name for internal consistency
#' @param fill_col Column name for heatmap fill values
#' @param fill_lab Legend label for fill scale
#' @return ggplot2 heatmap with significance annotation
#'
#' @details Design principles:
#' \itemize{
#'   \item Diverging color scale: Red (positive) to blue (negative) trends
#'   \item Symmetric limits: Equal emphasis on positive and negative changes
#'   \item Grey masking: Non-significant trends shown as neutral grey
#'   \item Significance overlay: Star symbols indicate confidence levels
#' }
create_trend_heatmap <- function(data, variable = NULL,
                                 fill_col = "sen_slope_plot",
                                 fill_lab = "Sen's slope") {
  # Calculate symmetric color scale limits based on data range
  trend_values <- data[[fill_col]]
  max_absolute_trend <- max(abs(trend_values), na.rm = TRUE)

  # Extract variable-specific legend label
  var_name <- as.character(unique(data$variable))
  legend_label <- var_labels[[var_name]]

  ggplot(data) +
    # Heatmap tiles colored by trend magnitude
    geom_tile(
      aes(
        x = factor(Station_Cluster), # X: Station clusters
        y = factor(classification), # Y: Water mass types
        fill = .data[[fill_col]] # Fill: Trend magnitude
      ),
      color = "white", # White tile borders for clarity
      linewidth = 0.5 # Thin border lines
    ) +

    # Significance symbols overlaid on tiles
    geom_text(
      aes(
        x = factor(Station_Cluster), # X: Match tile positions
        y = factor(classification), # Y: Match tile positions
        label = p_value_signif # Text: Significance symbols
      ),
      size = 7.5, # Large, readable font size
      color = "black", # White text for contrast
      fontface = "bold" # Bold for emphasis
    ) +

    # Water mass ordering (reversed for conventional oceanographic display)
    scale_y_discrete(limits = rev(Water_order)) + # Deep waters at bottom

    # Diverging color scale for trend direction and magnitude
    colorspace::scale_fill_continuous_divergingx(
      palette = "Geyser", # Professional diverging palette
      limits = c(-max_absolute_trend, max_absolute_trend), # Symmetric about zero
      name = fill_lab, # Legend title
      na.value = "grey75" # Non-significant trends in grey
    ) +

    # Axis labels and plot title
    labs(
      x = "Station", # Clear station identification
      y = "Water Mass", # Water mass classification
      title = legend_label # Variable-specific title
    ) +

    # Consistent theme for trend heatmaps
    base_theme
}

#' Assemble multi-panel trend heatmap figure with dynamic layout
#'
#' Creates comprehensive figure showing temporal trends across all measured
#' oceanographic variables in organized grid layout. Optimized for manuscript
#' Figure 4 with professional spacing and consistent labeling that reveals
#' coherent patterns of environmental change across the NOW transect.
#'
#' @param plot_function Function to generate individual trend heatmap panels
#' @param dataset Processed trend dataset containing Sen's slopes and significance
#' @param layout_config List containing variables, positions, and grid dimensions
#' @param fill_col Column name for heatmap coloring (trend magnitude)
#' @param fill_lab Legend label for color scale description
#' @return Combined ggplot2 object with multiple panels and individual legends
#'
#' @details Dynamic grid organization for temporal trend analysis:
#' \itemize{
#'   \item Row 1: Physical properties (CT, SA) establishing framework
#'   \item Row 2: Primary nutrients (Nitrate, Phosphate, Silicate, Oxygen) showing biogeochemical shifts
#'   \item Row 3: Derived tracers (N*) indicating water mass provenance changes
#'   \item Individual legends: Each panel retains scale for cross-variable comparison
#'   \item Flexible positioning: Accommodates different variable sets and layouts
#' }
#'
#' @examples
#' # Standard 7-variable layout for comprehensive trend analysis
#' layout_config <- list(
#'   variables = c("CT", "SA", "Nitrate", "Phosphate", "Silicate", "Oxconc", "Nstar"),
#'   positions = c(1, 2, 5, 6, 7, 8, 9),
#'   nrow = 3, ncol = 4
#' )
#' heatmap <- create_trend_grid(plot_function, trends, layout_config)
#'
create_trend_grid <- function(plot_function, dataset, layout_config,
                              fill_col = "sen_slope_plot",
                              fill_lab = "Sen's slope") {
  # Create blank plot for empty grid positions
  # Essential for maintaining visual balance in multi-panel scientific figures
  blank_plot <- ggplot() +
    theme_void()

  # Extract dynamic layout configuration parameters
  variables <- layout_config$variables # Oceanographic variables for trend analysis
  positions <- layout_config$positions # Strategic grid positions for each variable
  nrow <- layout_config$nrow # Number of rows for thematic organization
  ncol <- layout_config$ncol # Number of columns for balanced layout

  # Initialize grid structure with calculated dimensions
  total_slots <- nrow * ncol
  plots_list <- rep(list(blank_plot), total_slots) # Fill with blanks initially
  labels <- rep("", total_slots) # Empty labels for all positions

  # Generate trend heatmaps for each specified variable
  # Maintains consistent color scaling across panels for comparison
  for (i in 1:length(variables)) {
    pos <- positions[i]

    # Validate position within grid boundaries
    if (!is.na(pos) && pos <= total_slots) {
      # Filter dataset for current variable and generate heatmap
      # Each heatmap shows station×water_mass trends for one variable
      variable_data <- filter(dataset, variable == variables[i])

      plots_list[[pos]] <- plot_function(
        variable_data,
        variable = variables[i], # Pass variable name for plot customization
        fill_col = fill_col, # Column containing trend magnitudes
        fill_lab = fill_lab # Legend label for scientific interpretation
      )

      # Assign sequential panel labels for manuscript referencing
      labels[pos] <- LETTERS[i]
    }
  }

  # Assemble multi-panel trend figure with scientific formatting
  # Individual legends preserved to enable cross-variable trend comparison
  do.call(ggpubr::ggarrange, c(
    plots_list,
    list(
      labels = labels, # Panel labels (A, B, C, etc.) for manuscript
      nrow = nrow, # Grid dimensions from layout configuration
      ncol = ncol,
      align = "hv" # Align both horizontally and vertically
      # Note: Individual legends maintained for trend magnitude comparison
      # Unlike Figure 3, each variable needs its own scale for proper interpretation
    )
  ))
}

# ──────────────────────────────────────────────────────────────────────────────
# 4.4 GENERATE TEMPORAL TREND ANALYSIS FIGURE (MANUSCRIPT FIGURE 4)
# ──────────────────────────────────────────────────────────────────────────────

#' Create comprehensive temporal trend heatmap for decadal change analysis
#'
#' Multi-panel visualization revealing patterns of environmental change across
#' the North Water (NOW) transect from 2005-2023. Each panel shows station×water_mass
#' heatmaps with Sen's slope magnitudes and Mann-Kendall significance levels.
#'
#' Scientific objectives:
#' - Identify coherent spatial patterns of temporal change
#' - Distinguish variables experiencing strongest environmental shifts
#' - Reveal water-mass-specific responses to climate forcing
#' - Document evidence of Arctic Ocean transitions (Atlantification, Pacification)
#'
#' Analytical approach:
#' - Sen's slope: Robust trend estimator resistant to outliers
#' - Mann-Kendall test: Non-parametric significance assessment
#' - Station×water_mass framework: Links trends to oceanographic structure
#' - Color scaling: Red (warming/increasing), Blue (cooling/decreasing)

# Define layout configuration for comprehensive trend analysis
layout_config_trends <- list(
  # Variables ordered by oceanographic significance and interpretability
  variables = c("CT", "SA", "Nitrate", "Phosphate", "Silicate", "Nstar", "fPW"),

  # Strategic positioning for scientific narrative:
  # Row 1: Physical framework (T, S) + strategic spacing
  # Row 2: Complete nutrient suite (N, P, Si, O₂) showing biogeochemical shifts
  # Row 3: Water mass tracer (N*) + visual balance
  positions = c(1, 2, 4, 5, 6, 7, 8),

  # Grid dimensions optimized for trend comparison
  nrow = 3, # Three rows for thematic organization
  ncol = 3 # Four columns accommodate nutrient suite in single row
)

# Generate comprehensive temporal trend analysis figure
# This reveals decade-scale environmental changes across the NOW transect
temporal_trend_heatmap <- create_trend_grid(
  plot_function = create_trend_heatmap, # Function generating station×water_mass heatmaps
  dataset = temporal_trends_processed, # Complete trend dataset with Sen's slopes
  layout_config = layout_config_trends, # Layout configuration defined above
  fill_col = "sen_slope_plot", # Column containing trend magnitudes
  fill_lab = "Sen's slope (units/year)" # Scientific units for legend
)

# Display comprehensive trend analysis
# Expected to reveal coherent patterns of Arctic Ocean transition
temporal_trend_heatmap

# Export publication-ready trend figure with optimal specifications
# Dimensions accommodate multiple heatmaps with individual legends
# ggsave(
#   plot = temporal_trend_heatmap,
#   filename = paste0("wm_trends", file.format), # Descriptive manuscript filename
#   path = file.path, # Global output directory
#   width = 25, # Extra wide format for multi-panel heatmap layout
#   height = 10, # Sufficient height for three-row scientific presentation
#   dpi = 300 # Publication-quality resolution for journal submission
# )

# ──────────────────────────────────────────────────────────────────────────────
# 4.5 INDIVIDUAL TREND TIME SERIES ANALYSIS
# ──────────────────────────────────────────────────────────────────────────────

#' Prepare data for detailed time series visualization
#'
#' Combines trend statistics with raw observations to enable detailed
#' examination of individual station×water mass×variable combinations.
#' Useful for validating trend calculations and exploring temporal patterns.
#'
#' @details Data integration:
#' \itemize{
#'   \item Trend metadata: Sen's slope, p-values, significance symbols
#'   \item Raw observations: Individual measurements with decimal year
#'   \item Temporal coverage: Full time series from 2005-2023
#'   \item Quality filtering: Excludes mixed water classifications
#' }
#' individual_trend_data <- temporal_trends_processed %>%
#'   filter(!classification %in% "Mixed") %>% # Focus on classified water masses
#'   left_join(
#'     NOW %>%
#'       # Prepare raw observations for time series plotting
#'       select(
#'         Station_Cluster, eventDate, classification,
#'         CT, SA, Nitrate, Phosphate, Silicate, Oxconc, Nstar
#'       ) %>%
#'       pivot_longer(any_of(var_levels), names_to = "variable", values_to = "Value") %>%
#'       filter(!is.na(Value)) %>% # Remove missing observations
#'       mutate(
#'         # Convert dates to decimal years for trend line fitting
#'         dec_year = year(eventDate) + (yday(eventDate) - 1) /
#'           if_else(leap_year(eventDate), 366, 365)
#'       ),
#'     by = c("Station_Cluster", "classification", "variable") # Join on grouping variables
#'   )
#'
#' #' Generate individual trend plots for manuscript supplementary materials
#' plot_trend(individual_trend_data, 113, "PMW", "Nitrate")

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 5: PACIFIC WATER TRACER INTERCOMPARISON AND SENSITIVITY ANALYSIS
# ══════════════════════════════════════════════════════════════════════════════

# ──────────────────────────────────────────────────────────────────────────────
# 5.1 fPW METHOD SENSITIVITY ANALYSIS (continued)
# ──────────────────────────────────────────────────────────────────────────────

# ── 5.1.1 fPW estimates for each end-member / nitrogen-pool configuration ────
fPW_methods <- NOW %>%
  transmute(
    "Present study"                 = fpw(Phosphate, TIN,     atlantic_ref = "Jones2008",      pacific_ref = "Yamamoto.Kawai2008"),
    "Jones et al. (2008) - Nitrate" = fpw(Phosphate, Nitrate, atlantic_ref = "Jones2008",      pacific_ref = "Jones2008"),
    "Jones et al. (2008) - TIN"     = fpw(Phosphate, TIN,     atlantic_ref = "Jones2008",      pacific_ref = "Jones2008"),
    "Sutherland et al. (2009)"      = fpw(Phosphate, TIN,     atlantic_ref = "Sutherland2009", pacific_ref = "Sutherland2009"),
    "Sherwood et al. (2021)"        = fpw(Phosphate, TIN,     atlantic_ref = "Sherwood2021",   pacific_ref = "Sherwood2021")
  )

# Method names and line-broken labels (used by both table and plots)
methods       <- names(fPW_methods)
method_labels <- setNames(gsub(" \\(", "\n(", methods), methods)
n             <- length(methods)

# ── 5.1.2 Pairwise sensitivity statistics (unique pairs only) ────────────────
# combn() returns pairs (j, i) with j < i, i.e. the lower-triangle order
pairs <- t(combn(methods, 2))

fPW_sensitivity <- tibble(
  reference = pairs[, 1],   # x-axis of panel
  compared  = pairs[, 2]    # y-axis of panel
) %>%
  rowwise() %>%
  mutate(
    abs_diff      = list(abs(fPW_methods[[compared]] - fPW_methods[[reference]])),
    mean_abs_diff = mean(abs_diff, na.rm = TRUE),
    sd_abs_diff   = sd(abs_diff,   na.rm = TRUE),
    correlation   = cor(fPW_methods[[reference]], fPW_methods[[compared]],
                        use = "pairwise.complete.obs")
  ) %>%
  ungroup() %>%
  select(-abs_diff)

print(fPW_sensitivity)

# ── 5.1.3 Formatted table for the figure ─────────────────────────────────────
fPW_table <- fPW_sensitivity %>%
  transmute(
    "Reference method" = reference,
    "Compared method"  = compared,
    "Mean |ΔfPW|"      = sprintf("%.3f", mean_abs_diff),
    "SD |ΔfPW|"        = sprintf("%.3f", sd_abs_diff),
    "Pearson r"        = sprintf("%.3f", correlation)
  )

fPW_table_grob <- ggpubr::ggtexttable(
  fPW_table,
  rows  = NULL,
  theme = ggpubr::ttheme("light", base_size = 9)
)

# ── 5.1.4 Lower-triangle scatterplot matrix ──────────────────────────────────
plot_method_pair <- function(x_method, y_method) {
  ggplot(fPW_methods, aes(x = .data[[x_method]], y = .data[[y_method]])) +
    geom_point(size = 2, alpha = 0.5) +
    geom_abline(intercept = 0, slope = 1, col = "#CC0033", linewidth = 0.7) +
    labs(x = method_labels[[x_method]], y = method_labels[[y_method]]) +
    theme_bw(base_size = 10)
}

plot_matrix <- matrix(list(ggplot() + theme_void()), n, n)  # blank by default
for (i in 2:n) {
  for (j in 1:(i - 1)) {
    plot_matrix[[i, j]] <- plot_method_pair(methods[j], methods[i])
  }
}

# ggarrange fills row-wise -> transpose before flattening
fPW_method_comparison <- ggarrange(
  plotlist = as.vector(t(plot_matrix)),
  ncol = n, nrow = n,
  align = "hv"
)

# ── 5.1.5 Combined figure: matrix (A) + table (B) ────────────────────────────
fPW_method_figure <- ggarrange(
  fPW_method_comparison, fPW_table_grob,
  ncol   = 2,
  widths = c(1.6, 1),
  labels = c("A", "B")
)

# ggsave(
#   plot     = fPW_method_figure,
#   filename = paste0("fpw_method_comparison", file.format),
#   path     = file.path,
#   height   = 10, width = 20, dpi = 300
# )

# ──────────────────────────────────────────────────────────────────────────────
# 5.2 PACIFIC WATER TRACER CORRELATION MATRIX
# ──────────────────────────────────────────────────────────────────────────────

#' Evaluate relationships among Pacific Water tracers and assess method sensitivity
#' Supports tracer selection rationale and uncertainty quantification for manuscript

# Prepare tracer dataset with standardized variable names
pw_tracers <- NOW %>%
  select(Nstar, fPW, ANP, NO, PO, NO_PO, POs_star) %>%
  rename("N*" = "Nstar", "NO/PO" = "NO_PO", "POs*" = "POs_star")

# Custom correlation panel with color-coded background and significance stars
my_cor_panel <- function(data, mapping, ...) {
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  corr <- cor(x, y, use = "pairwise.complete.obs")

  # Color background based on correlation strength
  col <- colorRampPalette(c("#6D9EC1", "white", "#E46726"))(100)[
    findInterval(corr, seq(-1, 1, length.out = 100), all.inside = TRUE)
  ]

  # Significance testing for annotation
  p_val <- cor.test(x, y)$p.value
  stars <- symnum(p_val,
    corr = FALSE, na = FALSE,
    cutpoints = c(0, 0.001, 0.01, 0.05, 1),
    symbols = c("***", "**", "*", "")
  )

  ggplot() +
    theme_void() +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = col) +
    annotate("text",
      x = 0, y = 0,
      label = sprintf("%.3f%s", corr, stars),
      size = 5.5, fontface = "bold", color = "black"
    )
}

# Scatterplot panel with regression line
scatter_panel <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_point(alpha = 0.3, size = 1.5) +
    geom_smooth(method = "lm", col = "#CC0033", linewidth = 1.1, se = TRUE) +
    theme(panel.background = element_blank())
}

# Generate correlation matrix plot for tracer comparison
pw_correlation_matrix <- ggpairs(
  pw_tracers,
  lower = list(continuous = scatter_panel),
  diag = list(continuous = "blankDiag"),
  upper = list(continuous = my_cor_panel),
  progress = FALSE
) +
  theme(
    axis.text = element_text(size = 10, colour = "black"),
    axis.title = element_text(size = 15, colour = "black"),
    axis.line = element_line(colour = "black", linewidth = 0.75),
    axis.ticks = element_line(colour = "black", linewidth = 0.75),
    strip.text = element_text(size = 17.5, face = "bold")
  )

pw_correlation_matrix

# Save high-resolution figure for supplementary materials
# ggsave(
#   plot = pw_correlation_matrix,
#   filename = paste0("pw_correlation_matrix", file.format),
#   path = file.path,
#   height = 10, width = 15, dpi = 300
# )

# ──────────────────────────────────────────────────────────────────────────────
# 5.4 PACIFIC WATER FRACTION DISTRIBUTION BY WATER MASS
# ──────────────────────────────────────────────────────────────────────────────

# Analysis of Pacific Water fraction (fPW) distribution across water masses
# Higher fPW values indicate stronger Pacific Ocean influence
fpw_boxplot <- ggplot(data = NOW, aes(y = classification, x = fPW, fill = classification)) +
  geom_boxplot(
    outlier.shape = 21, outlier.size = 2, outlier.fill = "white", outlier.stroke = 1.5,
    linewidth = 1.5
  ) +
  # geom_jitter(width = 0, height = 0.15, alpha = 0.3, size = 0.5) +
  scale_fill_manual(values = wm_palette, name = "Water Mass") +
  scale_y_discrete(limits = rev(Water_order)) +
  labs(
    y = "Water Mass",
    x = "fPW",
    # title = "Distribution of Pacific Water Influence Across Water Masses"
  ) +
  plot_theme +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 2),
    axis.line = element_blank(),
    legend.position = "none"
  )

print(fpw_boxplot)

# Save high-resolution figure for supplementary materials
ggsave(
  plot = fpw_boxplot,
  filename = paste0("fpw_boxplot", file.format),
  path = file.path,
  height = 7.5, width = 10, dpi = 300
)

# ──────────────────────────────────────────────────────────────────────────────
# 5.5 SPATIAL DISTRIBUTION OF PACIFIC WATER FRACTION
# ──────────────────────────────────────────────────────────────────────────────

# Spatial distribution of median fPW values across the NOW transect
# Reveals geographic patterns of Pacific vs Atlantic water mass influence
# Calculate spatial distribution with additional statistics
spatial_fpw <- NOW %>%
  filter(!is.na(fPW)) %>% # Remove missing fPW values
  group_by(Station_Cluster, classification) %>%
  summarise(
    median_fPW = median(fPW, na.rm = TRUE), # Central tendency
    mean_fPW = mean(fPW, na.rm = TRUE), # Alternative measure
    iqr_fPW = IQR(fPW, na.rm = TRUE), # Variability
    n = n(), # Sample size per cell
    .groups = "drop"
  )

# Create annotated heatmap
fpw_spatial_distrib <- spatial_fpw %>%
  ggplot(aes(x = factor(Station_Cluster), y = classification, fill = median_fPW)) +
  geom_tile(color = "white", size = 0.75, alpha = 0.9) +
  geom_text(aes(label = round(median_fPW, 2)), size = 3, fontface = "bold") +

  # Color scale for Pacific water fraction
  scale_fill_gradient2(
    low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
    midpoint = 0.5, name = "fPW",
    limits = c(0, 1)
  ) +

  # Labels and annotations
  labs(
    # title = "Spatial Distribution of Pacific Water Fraction Along NOW Transect",
    y = "Water Mass",
    x = "Station"
  ) +
  scale_y_discrete(limits = rev(Water_order)) +
  # Theme adjustments
  plot_theme +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 2),
    axis.line = element_blank()
  )

print(fpw_spatial_distrib)

# Save high-resolution figure for supplementary materials
# ggsave(
#   plot = fpw_spatial_distrib,
#   filename = paste0("fpw_spatial_distribution", file.format),
#   path = file.path,
#   height = 5, width = 10, dpi = 300
# )

# ──────────────────────────────────────────────────────────────────────────────
# 5.6 CROSS-SECTIONAL WATER MASS DISTRIBUTION WITH PACIFIC INFLUENCE
# ──────────────────────────────────────────────────────────────────────────────

#' Generate multi-year composite section showing water mass distribution
#' with Pacific Water fraction isolines to illustrate Arctic-Pacific connectivity
wm_section <- plot_watermass_section(
  NOW,
  contour_variable = "fPW", # Pacific Water fraction isolines
  contour_breaks = c(0.5, 0.7), # Moderate (50%) and strong (70%) Pacific influence
  ngrid = 75 # High-resolution spatial interpolation
)

wm_section

# Save high-resolution figure for manuscript
# ggsave(
#   plot = wm_section,
#   filename = paste0("wm_section", file.format),
#   path = file.path,
#   width = 20,
#   height = 10,
#   dpi = 300
# )

# ──────────────────────────────────────────────────────────────────────────────
# 5.7 CROSS-SECTIONAL fPW DISTRIBUTION
# ──────────────────────────────────────────────────────────────────────────────
fpw_section <- plot_section(NOW, param = "fPW", contour_var = "self")

# --- Display t   he plot
fpw_section

# Save high-resolution figure for manuscript
# ggsave(
#   plot = fpw_section,
#   filename = paste0("fpw_section", file.format),
#   path = file.path,
#   width = 22.5,
#   height = 12.5,
#   dpi = 300
# )

# ══════════════════════════════════════════════════════════════════════════════
# 5.8 PACIFIC WATER REGIME CLASSIFICATION
# ══════════════════════════════════════════════════════════════════════════════
pw_regime_df <- NOW %>% filter(Depth <= 200, !is.na(fPW))

# pw_regime_df %>%
#   count(classification, sort = T) %>%
#   mutate(prop = round(n / sum(n) * 100, 0))

pw_regime_results <- pw_regime_calc(
  pw_regime_df,
  min_stations        = 5,
  fPW_cutoff          = 0.50,
  front_detection_thr = 0.50,
  western_limit       = -74,
  surface_only        = FALSE
)

# >>> SAMPLING-TIME SUMMARY (for figure captions or legends)
sampling_meta <- NOW %>%
  group_by(year) %>%
  summarise(
    mid_date = as.Date((min(DayOfYear) + max(DayOfYear)) / 2 - 1,
      origin = paste0(year[1], "-01-01")
    ),
    .groups = "drop"
  ) %>%
  mutate(
    period = case_when(
      as.integer(format(mid_date, "%d")) <= 10 ~ "Early",
      as.integer(format(mid_date, "%d")) <= 20 ~ "Mid",
      TRUE ~ "Late"
    ),
    month_period = paste(period, month.abb[as.integer(format(mid_date, "%m"))])
  ) %>%
  select(year, month_period)


# ──────────────────────────────────────────────────────────────────────────────
# REGIME CLASSIFICATION ANALYSIS
# ──────────────────────────────────────────────────────────────────────────────
# Apply 6-tier classification system:
# 1. Atlantic Dominance         - Complete Pacific collapse (< 5% influence)
# 2. Marginal Pacific Inflow    - Minimal Pacific presence (5-25% influence)
# 3. Mixed Dominance            - Transitional state (25-40% influence)
# 4. Western Accumulation       - High volume but spatially confined (≥40%, West of -74°W)
# 5. Pacific Expansion          - Moderate eastward extension (40-60% influence)
# 6. Pacific Saturation         - Extreme Pacific dominance (≥60% influence)

pw_regime_table <- pw_regime_results %>%
  left_join(sampling_meta, by = "year") %>%
  select(year, month_period, QC_Flag, Regime, Dominance_Index, Front_Lon) %>%
  mutate(
    Dominance_Index = if_else(QC_Flag == "Exclude", NA, round(Dominance_Index, 2)),
    Front_Lon = round(Front_Lon, 1)
  )

# ── Single-table export ──────────────────────────────────────────────────────
# writexl::write_xlsx(
#   x = pw_regime_table, # Distance summary table
#   path = file.path(
#     file.path, # Global output directory
#     paste0("pw_regime_table", ".xlsx")
#   ) # Descriptive filename
# )

# ──────────────────────────────────────────────────────────────────────────────
# SUPPLEMENTARY FIGURE — Si/P RATIO (SPATIAL & TEMPORAL SYNTHESIS)
# ──────────────────────────────────────────────────────────────────────────────

#' Panel A — Temporal trends of the Si/P ratio
#'
#' Heatmap of Sen's slopes for the Si/P ratio across stations and water masses.
#' Data are filtered to retain the single derived tracer of interest (SiP)
#' before being passed to the generic trend-heatmap constructor.
SiP_Sen <-
  temporal_trends_processed %>%
  filter(variable == "SiP") %>% # Keep Si/P tracer only
  create_trend_heatmap(
    fill_col = "sen_slope_plot", # Fill: Sen's slope magnitude
    fill_lab = "Sen's slope" # Legend title
  )

#' Panel B — Spatial distribution of the Si/P ratio
#'
#' Line-range plot summarizing the median Si/P ratio (with IQR) by water mass
#' along the NOW transect. Longitude axis enabled to anchor the spatial reading.
SiP_wm <- plot_linerange(
  wm_stat %>% filter(variable %in% "SiP"), # Subset Si/P statistics
  show_lon_axis = TRUE # Display longitude reference
) + theme(legend.position = "right")

# ── Compose the two-panel supplementary figure ───────────────────────────────
#' Combine temporal (A) and spatial (B) views into one publication-ready figure.
SiP <- ggpubr::ggarrange(
  SiP_wm, SiP_Sen, # Two main visualization panels
  # get_legend,                                                          # (optional shared legend)
  labels = LETTERS, # Scientific panel labels (A, B)
  nrow = 1, # Single-row layout
  ncol = 2, # Side-by-side columns
  widths = c(1, 1), # Equal panel widths
  heights = c(1, 1), # Equal panel heights
  align = "v" # Vertical alignment for clean look
)

# ── Export ────────────────────────────────────────────────────────────────────
#' Write the supplementary Si/P figure to disk at publication quality.
# ggsave(
#   plot     = SiP,
#   filename = paste0("SiP", file.format), # Descriptive filename
#   path     = file.path, # Global output directory
#   width    = 20, # Wide format for side-by-side panels
#   height   = 5, # Appropriate height for detail
#   dpi      = 300 # Publication-quality resolution
# )

SPMW <- NOW %>%
  filter(classification == "SPMW",
         !is.na(SiP), !is.na(fPW), is.finite(SiP))

n_obs <- nrow(SPMW)

SPMW_SiP_fPW <- 
  ggplot(SPMW, aes(x = SiP, y = fPW)) +
  geom_point(shape = 21,fill = "#B276B2FF", size = 2.5, stroke = 1) +
  geom_smooth(method = "lm", formula = y ~ x, colour = "black", linewidth = 0.8) +
  stat_cor(method = "spearman", cor.coef.name = "rho",
           aes(label = paste(after_stat(r.label), after_stat(p.label),
                             paste0("italic(n)~`=`~", n_obs),
                             sep = "~`,`~")),
           label.x.npc = 0.70,                       # 0 = left, 1 = right
           label.y.npc = 0.98,
           r.digits = 2, p.accuracy = 0.001, size = 3.5) +
  # coord_cartesian(xlim = c(0, 25), ylim = c(0, 1)) +
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.15))) +
  scale_x_continuous(expand = expansion(mult = c(0.15, 0.15))) +
  labs(x = "Si:P", y = expression(italic(f)[PW])) +
  plot_theme +
  theme(
    legend.position = "none",
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 2),
    strip.background = element_rect(fill = "grey50"),
    strip.text       = element_text(face = "bold"),
    axis.line        = element_blank() # Clean borders instead of axis lines
  )

# ── Export ────────────────────────────────────────────────────────────────────
#' Write the supplementary Si/P figure to disk at publication quality.
ggsave(
  plot     = SPMW_SiP_fPW,
  filename = paste0("SPMW_SiP_fPW", file.format), # Descriptive filename
  path     = file.path, # Global output directory
  width    = 10, # Wide format for side-by-side panels
  height   = 5, # Appropriate height for detail
  dpi      = 300 # Publication-quality resolution
)
