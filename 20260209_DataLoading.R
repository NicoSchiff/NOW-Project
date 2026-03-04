# ══════════════════════════════════════════════════════════════════════════════
# DATA LOADING AND PREPROCESSING
# ══════════════════════════════════════════════════════════════════════════════

#' Load and combine multi-year nutrient datasets from ArcticNet expeditions
#'
#' This function processes and standardizes nutrient datasets from multiple years
#' ensuring consistent temporal and spatial referencing across all observations.
#'
#' @param parameters character vector of biogeochemical parameters to extract
#' @return data.frame with standardized multi-year nutrient observations
#'
#' @details
#' Data sources include:
#' - R_Database_dec2021.txt: Historical data (2005-2021) 
#' - R_Data_2022.txt: 2022 expedition data
#' - R_Data_2023.txt: 2023 expedition data
#'
#' Standardization procedures:
#' - Forward filling of station metadata
#' - Longitude convention conversion (Western positive → Eastern negative)
#' - Temporal component extraction and day-of-year calculation
#' - Duplicate removal and data validation
load_data <- function(parameters = parameter) {
  
  # Define source files containing hydrographic and nutrient observations
  files <- c(
    "DataBase_Nut/Data/R_Database_dec2021.txt",  # Historical data through 2021
    "DataBase_Nut/Data/R_Data_2022.txt",         # 2022 expedition data
    "DataBase_Nut/Data/R_Data_2023.txt"          # 2023 expedition data
  )
  
  # Process each dataset independently to maintain data integrity
  processed_datasets <- map(files, function(file_path) {
    cat("Processing:", basename(file_path), "\n")
    
    # Load raw data with empty strings treated as NA
    raw_data <- read.delim(file_path, na.strings = "")
    
    # Apply consistent data transformations
    processed_data <- raw_data %>%
      # Forward fill station metadata to handle missing values in CTD profiles
      fill(Cruise:Bot.Depth) %>%
      mutate(
        # Parse and extract temporal components for time series analysis
        date_parsed = as.POSIXct(eventDate, format = "%Y-%m-%d"),
        year = year(date_parsed),
        month = month(date_parsed),
        day = day(date_parsed),
        DayOfYear = as.numeric(strftime(make_date(year, month, day), format = "%j")),
        
        # Convert longitude from Western (positive) to Eastern (negative) convention
        # This ensures consistency with oceanographic mapping standards
        decimalLongitude = decimalLongitude - 360
      )
    
    # Standardize column selection across all datasets
    selected_data <- processed_data %>%
      select(
        # Core expedition and station identification
        Cruise, Station, Type, eventDate, year, month, day, DayOfYear,
        
        # Geographic positioning (WGS84 coordinates)
        decimalLongitude, decimalLatitude, Depth,
        
        # Water mass stratification information (if available in dataset)
        any_of("Strate"),
        
        # Target biogeochemical parameters for NOW transect analysis
        any_of(parameters)
      ) %>%
      arrange(year, month, day, Cruise, Station) %>%
      distinct()  # Remove potential duplicates from overlapping datasets
    
    return(selected_data)
  })
  
  # Combine all processed datasets into single comprehensive dataframe
  final_data <- bind_rows(processed_datasets)
  
  # Generate data loading summary for validation
  cat("Data loading complete!\n")
  cat("Total observations:", nrow(final_data), "\n")
  cat("Years covered:", paste(range(final_data$year, na.rm = TRUE), collapse = " - "), "\n")
  cat("Geographic extent:\n")
  cat("  Longitude:", paste(round(range(final_data$decimalLongitude, na.rm = TRUE), 2), collapse = " to "), "\n")
  cat("  Latitude:", paste(round(range(final_data$decimalLatitude, na.rm = TRUE), 3), collapse = " to "), "\n")
  
  return(final_data)
}
