# ──────────────────────────────────────────────────────────────────────────────
# FONCTION MULTI-ANNÉES pour sections de masses d'eau
# ──────────────────────────────────────────────────────────────────────────────

plot_watermass_section <- function(
    data, 
    year_select = NULL,
    ngrid = 100,
    lambda = 0.005,
    threshold = 0.6,
    comparable_ratio = 1.3,
    add_obs_points = TRUE,
    add_contours = TRUE,           
    contour_variable = "fPW",      
    contour_breaks = c(0.3, 0.5, 0.7),
    bathy_file = "NOW_bathymetry.csv",
    endmembers = wm_ref
) {
  library(tidyverse)
  library(fields)
  library(ggplot2)
  library(metR)
  library(pals)
  library(MBA)
  library(reshape2)  # ← IMPORTANT pour melt()
  
  # ────────────────────────────────────────────────────────────────────────────
  # 0) Bathymétrie
  # ────────────────────────────────────────────────────────────────────────────
  
  bathy_data <- read.csv("/Users/nicolas/Library/CloudStorage/OneDrive-UniversitéLaval/DataBase_Nut/Data/NOW_bathymetry.csv")
  
  # ────────────────────────────────────────────────────────────────────────────
  # 1) Préparation des données
  # ────────────────────────────────────────────────────────────────────────────
  
  if (!is.null(year_select)) {
    data <- data %>% filter(year %in% year_select)
  }
  
  clean_data <- data %>% 
    filter(!is.na(sigma0), !is.na(spicity), !is.na(decimalLongitude), !is.na(Depth)) %>%
    arrange(year, decimalLongitude, Depth)
  
  # ────────────────────────────────────────────────────────────────────────────
  # 2) CONTOURS - TON APPROCHE QUI MARCHE ! 🎯
  # ────────────────────────────────────────────────────────────────────────────
  
  interp_contour_grid <- NULL
  contour_linetypes <- NULL
  
  if (add_contours && contour_variable %in% names(clean_data)) {
    
    # ──── TON CODE ORIGINAL ──── 
    interp_contour <- clean_data %>%
      select(decimalLongitude, Depth, year, contains(contour_variable)) %>%
      drop_na(contains(contour_variable)) %>%
      group_by(year) %>%
      group_split()
    
    interp_contour_raster <- lapply(interp_contour, function(df) {
      if (nrow(df) < 5) {
        return(NULL)
      }
      
      tryCatch({
        mba <- mba.surf(df[, c("decimalLongitude", "Depth", contour_variable)], 
                        no.X = 100, no.Y = 100, extend = TRUE)
        zmat <- mba$xyz.est$z
        dimnames(zmat) <- list(mba$xyz.est$x, mba$xyz.est$y)
        
        reshape2::melt(zmat, varnames = c("decimalLongitude", "Depth"), 
                       value.name = contour_variable) %>%
          mutate(year = df$year[1]) %>%
          filter(!is.na(.data[[contour_variable]]))
      }, error = function(e) {
        cat("Error for year", unique(df$year), ":", e$message, "\n")
        return(NULL)
      })
    })
    
    # Combiner les résultats
    interp_contour_grid <- bind_rows(compact(interp_contour_raster))
    
    if (!is.null(interp_contour_grid) && nrow(interp_contour_grid) > 0) {
      cat("✅ Contour data created:", nrow(interp_contour_grid), "points\n")
      
      # Styles de lignes
      contour_linetypes <- setNames(
        rep(c("dotted", "solid", "dashed"), length.out = length(contour_breaks)),
        as.character(contour_breaks)
      )
    } else {
      add_contours <- FALSE
      cat("❌ No contour data created\n")
    }
  }
  
  # ────────────────────────────────────────────────────────────────────────────
  # 3) Classification des masses d'eau (ton code existant...)
  # ────────────────────────────────────────────────────────────────────────────
  
  get_seafloor_depth <- function(longitude) {
    approx(x = bathy_data$decimalLongitude, y = bathy_data$Depth, 
           xout = longitude, rule = 2)$y
  }
  
  interpolate_with_bathy <- function(data_subset, variable) {
    tps_fit <- Tps(
      x = data_subset[, c("decimalLongitude", "Depth")],
      Y = data_subset[[variable]],
      lambda = lambda
    )
    
    lon_range <- range(data_subset$decimalLongitude)
    lon_seq <- seq(lon_range[1], lon_range[2], length.out = ngrid)
    depth_seq <- seq(0, max(bathy_data$Depth, na.rm = TRUE), length.out = ngrid)
    pred_grid <- expand_grid(decimalLongitude = lon_seq, Depth = depth_seq)
    
    pred_values <- predict(tps_fit, pred_grid[, c("decimalLongitude", "Depth")])
    
    result <- pred_grid %>%
      mutate(
        Value = as.vector(pred_values),
        seafloor_depth = get_seafloor_depth(decimalLongitude)
      ) %>%
      filter(Depth <= seafloor_depth, !is.na(Value))
    
    return(result)
  }
  
  # Classification par année...
  all_classified <- list()
  
  for (selected_year in sort(unique(clean_data$year))) {
    cat("\nProcessing year:", selected_year, "\n")
    
    data_year <- clean_data %>% filter(year == selected_year)
    
    sigma0_grid <- interpolate_with_bathy(data_year, "sigma0") %>% 
      rename(sigma0 = Value) %>% select(-seafloor_depth)
    
    spicity_grid <- interpolate_with_bathy(data_year, "spicity") %>% 
      rename(spicity = Value) %>% select(-seafloor_depth)
    
    combined_grid <- sigma0_grid %>% 
      full_join(spicity_grid, by = c("decimalLongitude", "Depth")) %>% 
      filter(!is.na(sigma0) & !is.na(spicity)) %>%
      mutate(year = selected_year)
    
    classified_grid <- classify_watermass(
      data = combined_grid,
      endmembers = wm_ref,
      threshold = threshold,
      comparable_ratio = comparable_ratio,
      return_long = FALSE
    )
    
    all_classified[[as.character(selected_year)]] <- classified_grid
  }
  
  final_classified <- bind_rows(all_classified)
  
  # ────────────────────────────────────────────────────────────────────────────
  # 4) Visualisation
  # ────────────────────────────────────────────────────────────────────────────
  
  watermass_colors <- c(
    "PMW"  = "#5DA5DAFF",  
    "UHW"  = "#FAA43AFF",  
    # "CBAW" = "#F15854FF",  
    "BBPW" = "#60BD68FF",  
    "SPMW" = "#B276B2FF", 
    "Mixed" = "grey25"
  )
  
  section.WaterMass_plot <- ggplot()
  
  for (mass in names(watermass_colors)) {
    subset_mass <- subset(final_classified, classification == mass)
    if(nrow(subset_mass) > 0) {
      section.WaterMass_plot <- section.WaterMass_plot +
        ggnewscale::new_scale_fill() +
        geom_raster(
          data = subset_mass,
          aes(x = decimalLongitude, y = Depth, fill = distance)
        ) +
        scale_fill_gradient(
          low = watermass_colors[mass],
          high = "white",
          limits = c(0, threshold),
          name = paste(mass, "dist."),
          na.value = "grey50",
          guide = guide_coloursteps(
            direction = "horizontal",
            title.position = "top",
            title.hjust = 0.5,
            barwidth = unit(3, "cm"),
            barheight = unit(0.5, "cm")
          )
        )
    }
  }
  
  # ──── AJOUT DES CONTOURS - TON APPROCHE ! ────
  if (add_contours && !is.null(interp_contour_grid)) {
    section.WaterMass_plot <- section.WaterMass_plot +
      geom_contour(
        data = interp_contour_grid,
        aes(
          x = decimalLongitude,
          y = Depth,
          z = .data[[contour_variable]], # ← TON SYNTAXE QUI MARCHE
          linetype = factor(after_stat(level))  # ← UPDATED SYNTAX
        ),
        breaks = contour_breaks,
        color = "black"
      ) +
      scale_linetype_manual(
        values = contour_linetypes,
        name = paste0(contour_variable, " level")
      ) +
      metR::geom_text_contour(
        data = interp_contour_grid,
        aes(x = decimalLongitude, y = Depth, z = .data[[contour_variable]]), # ← TON SYNTAXE
        breaks = contour_breaks,
        stroke.colour = "white", 
        stroke = 0.25, 
        color = "black", 
        skip = 0
      )
  }
  
  # ──── FINALISATION ────
  p_final <- section.WaterMass_plot +
    geom_polygon(
      data = bathy_data %>%
        bind_rows(
          tibble(decimalLongitude = max(bathy_data$decimalLongitude),
                 Depth = max(bathy_data$Depth, na.rm = TRUE)),
          tibble(decimalLongitude = min(bathy_data$decimalLongitude),
                 Depth = max(bathy_data$Depth, na.rm = TRUE))
        ),
      aes(x = decimalLongitude, y = Depth), 
      fill = "grey20", color = "black", size = 0.8
    ) +
    scale_x_continuous(
      expand = c(0, 0),
      breaks = seq(-78, -72, 2),
      labels = function(x) paste0(abs(x), "°W")
    ) +
    scale_y_reverse(expand = c(0, 0)) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(size = 11, colour = "black"),
      title = element_text(size = 12, face = "bold"),
      legend.position = "right",
      strip.text = element_text(size = 12, face = "bold")
    )
  
  if (add_obs_points) {
    p_final <- p_final +
      geom_point(
        data = clean_data, 
        aes(x = decimalLongitude, y = Depth),
        size = 0.6, color = "white", alpha = 0.9
      )
  }
  
  if (length(unique(final_classified$year)) > 1) {
    p_final <- p_final + facet_wrap(~ year)
  }
  
  p_final <- p_final +
    labs(
      title = "Water Mass Classification - NOW Transect",
      x = "Longitude", 
      y = "Depth (m)"
    )
  
  return(p_final)
}

# Test
# p_with_fPW <- plot_watermass_section(NOW, contour_variable = "fPW")
