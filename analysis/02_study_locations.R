library(tidyverse)
library(lubridate)
library(janitor)

# Function to parse coordinate strings with various formats
parse_coords <- function(coord_string) {
  if (is.na(coord_string)) return(tibble(latitude = NA_real_, longitude = NA_real_))

  # First try to extract decimal coordinates in parentheses
  decimal_coords <- str_extract_all(coord_string, "-?\\d+\\.\\d+")[[1]]
  if (length(decimal_coords) >= 2) {
    # Convert to numeric and handle direction (N/S/E/W)
    lat <- as.numeric(decimal_coords[1])
    lon <- as.numeric(decimal_coords[2])

    # Check for S/W directions which make the coordinate negative
    if (str_detect(coord_string, "[sS]")) lat <- -abs(lat)
    if (str_detect(coord_string, "[wW]")) lon <- -abs(lon)

    return(tibble(latitude = lat, longitude = lon))
  }

  # If no decimal coordinates found, try to parse DMS format
  dms_pattern <- "(\\d+)°(\\d+)?'?(\\d+)?\"?([NSEWnsew])?"
  dms_matches <- str_match_all(coord_string, dms_pattern)[[1]]

  if (nrow(dms_matches) >= 2) {
    dms_to_decimal <- function(dms) {
      deg <- as.numeric(dms[2])
      min <- ifelse(is.na(dms[3]), 0, as.numeric(dms[3])/60)
      sec <- ifelse(is.na(dms[4]), 0, as.numeric(dms[4])/3600)
      dir <- toupper(ifelse(is.na(dms[5]), "", dms[5]))

      decimal <- deg + min + sec
      if (dir %in% c("S", "W")) decimal <- -decimal
      return(decimal)
    }

    lat <- dms_to_decimal(dms_matches[1,])
    lon <- dms_to_decimal(dms_matches[2,])
    return(tibble(latitude = lat, longitude = lon))
  }

  return(tibble(latitude = NA_real_, longitude = NA_real_))
}

# Function to parse date ranges
parse_date_range <- function(date_str) {
  if (is.na(date_str)) return(list(date_start = NA, date_end = NA))

  # Extract all years (4-digit numbers)
  years <- as.numeric(str_extract_all(date_str, "\\b\\d{4}\\b")[[1]])

  if (length(years) == 0) {
    return(list(date_start = NA, date_end = NA))
  } else if (length(years) == 1) {
    return(list(
      date_start = as.Date(paste0(years[1], "-01-01")),
      date_end = as.Date(paste0(years[1], "-12-31"))
    ))
  } else {
    return(list(
      date_start = as.Date(paste0(min(years), "-01-01")),
      date_end = as.Date(paste0(max(years), "-12-31"))
    ))
  }
}

# Read the data with proper encoding
raw_data <- read_csv("analysis/data-raw/study_locations.csv",
                    locale = locale(encoding = "latin1"),
                    trim_ws = TRUE) %>%
  # Clean column names and handle the combined lat/long column
  rename_with(~ "lat_long", matches("Latitude|Longitude")) %>%
  rename_with(tolower) %>%
  janitor::clean_names() %>%
  mutate(across(where(is.character), ~ str_trim(.x, "both")))

# Process the data
processed_data <- raw_data %>%
  # First clean up the lat_long column
  mutate(lat_long = str_replace_all(lat_long, "\\s*;\\s*", "\\n")) %>%
  # Split into multiple rows if there are multiple locations
  separate_rows(lat_long, sep = "\\n\\s*\\n") %>%
  filter(!is.na(lat_long) & lat_long != "") %>%
  # Extract coordinates and dates
  rowwise() %>%
  mutate(
    coords = list(parse_coords(lat_long)),
    date_info = list(parse_date_range(year_s))
  ) %>%
  unnest_wider(c(coords, date_info)) %>%
  # Clean up column order and names
  select(
    study,
    location = lat_long,
    latitude,
    longitude,
    date_start,
    date_end,
    mean_age_yrs = mean_age_yrs,
    age_group = children_adults_both,
    everything()
  )

# View the processed data
head(processed_data)

# manual cleaning
processed_data$latitude[18] <- 14.2
processed_data$longitude[18] <- 104.1
processed_data$latitude[20] <- 5.56
processed_data$longitude[20] <- -0.19
processed_data$latitude[12] <- 23.7
processed_data$longitude[9] <- 43.72

# TODO: fix the years properly

# get the prevalence
# sort out middate range for getting prevalence
processed_data <- processed_data %>%
  mutate(mid_date = date_start + (date_end - date_start)/2)

# Lastly add in MAP estimated malaria prevalence for these regions
prevs <- processed_data %>%
  split(processed_data$location) %>%
  map(function(x){
    rast <- malariaAtlas::getRaster("Malaria__202406_Global_Pf_Parasite_Rate",
                                    year = unique(lubridate::year(x$mid_date[1])),
                                    extent = matrix(c(min(x$longitude)-0.2, min(x$latitude)-0.2, max(x$longitude)+0.2, max(x$latitude)+0.2), nrow = 2))
    points <- x %>% select(latitude, longitude)
    pointst <- terra::vect(points, geom = c("longitude", "latitude"), crs = "EPSG:4326")

    # Extract raster values at these points
    extracted_values <- terra::extract(rast, pointst)
    if(is.na(extracted_values[1,2])) {
      x$prev <- mean(rast[,][,1], na.rm = TRUE)
    } else {
    x$prev <- as.numeric(extracted_values[,2])
    }
    return(x)
  })

final <- do.call(rbind, prevs)


# Create output directory if it doesn't exist
if (!dir.exists("analysis/data-derived")) {
  dir.create("analysis/data-derived", recursive = TRUE)
}

# Save the processed data
write_csv(final, "analysis/data-derived/study_locations_processed.csv")
