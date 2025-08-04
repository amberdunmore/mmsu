# Extract prevalence data for missing studies: Ashley 2014 and Chang 2016

library(malariaAtlas)
library(terra)
library(tidyverse)
library(lubridate)

# Create missing studies data based on your table information
missing_studies <- data.frame(
  study = c(
    # Ashley 2014 - 7 populations across different countries
    rep("Ashley 2014", 7),
    # Chang 2016 - Uganda location
    rep("Chang 2016", 2)  # You have 2 Chang 2016 entries in malaria_data
  ),
  location = c(
    # Ashley 2014 locations
    "Cambodia (Pailin)", "Cambodia (Preah Vihear)", "Cambodia (Ratanakiri)", "Cambodia (Pursat)",
    "Thailand (Mae Sot)", "Thailand (Srisaket)", "Thailand (Ranong)",
    # Chang 2016 - Uganda (general location since exact coordinates not provided)
    "Uganda (general)", "Uganda (general)"
  ),
  country = c(
    rep("Cambodia", 4), rep("Thailand", 3), rep("Uganda", 2)
  ),
  # Approximate coordinates based on the locations mentioned
  latitude = c(
    # Cambodia sites (from your Ashley table)
    12.85, 13.80, 13.73, 12.53,  # Cambodia coordinates from table
    # Thailand sites
    16.72, 15.12, 9.97,          # Thailand coordinates from table
    # Uganda - use country centroid since exact location not provided
    1.0, 1.0
  ),
  longitude = c(
    # Cambodia sites
    102.60, 104.97, 107.00, 103.92,  # Cambodia coordinates
    # Thailand sites
    98.57, 104.32, 98.63,            # Thailand coordinates
    # Uganda
    32.0, 32.0
  ),
  # Study dates
  date_start = c(
    rep(as.Date("2011-01-01"), 7),  # Ashley 2014 - approximate study period
    rep(as.Date("2013-01-01"), 2)   # Chang 2016 - approximate study period
  ),
  date_end = c(
    rep(as.Date("2013-12-31"), 7),  # Ashley 2014
    rep(as.Date("2015-12-31"), 2)   # Chang 2016
  ),
  stringsAsFactors = FALSE
)

print("Missing studies to extract prevalence for:")
print(missing_studies)

# Calculate mid-date for MAP extraction
missing_studies <- missing_studies %>%
  mutate(mid_date = date_start + (date_end - date_start)/2)

# Extract prevalence using the same method as your professor
cat("\n=== EXTRACTING PREVALENCE DATA ===\n")

extract_prevalence_for_studies <- function(study_data) {
  results_list <- list()

  for(i in 1:nrow(study_data)) {
    row <- study_data[i, ]
    cat("Extracting for:", row$study, "-", row$location, "\n")

    tryCatch({
      # Get MAP raster for the study year
      year <- year(row$mid_date)
      cat("  Year:", year, "\n")

      # Get raster with buffer around coordinates
      rast <- malariaAtlas::getRaster(
        "Malaria__202406_Global_Pf_Parasite_Rate",
        year = year,
        extent = matrix(c(
          row$longitude - 0.2, row$latitude - 0.2,
          row$longitude + 0.2, row$latitude + 0.2
        ), nrow = 2)
      )

      # Create point for extraction
      point <- terra::vect(
        data.frame(lon = row$longitude, lat = row$latitude),
        geom = c("lon", "lat"),
        crs = "EPSG:4326"
      )

      # Extract prevalence
      extracted_values <- terra::extract(rast, point)

      if(is.na(extracted_values[1,2])) {
        # If point extraction fails, use regional average
        cat("  Point extraction failed, using regional average\n")
        prev_value <- mean(terra::values(rast), na.rm = TRUE)
      } else {
        prev_value <- as.numeric(extracted_values[,2])
      }

      cat("  Extracted prevalence:", prev_value, "\n")

      # Store result
      results_list[[i]] <- row %>%
        mutate(prev = prev_value)

    }, error = function(e) {
      cat("  Error extracting for", row$study, "-", row$location, ":", e$message, "\n")
      # Use a default value or skip
      results_list[[i]] <- row %>%
        mutate(prev = NA)
    })
  }

  return(do.call(rbind, results_list))
}

# Extract prevalence data
new_prevalence_data <- extract_prevalence_for_studies(missing_studies)

# Show results
cat("\n=== EXTRACTION RESULTS ===\n")
print(new_prevalence_data %>%
        select(study, location, latitude, longitude, prev) %>%
        mutate(prev_pct = prev * 100))

# Read existing processed data
existing_data <- read_csv("analysis/data-derived/study_locations_processed.csv")

# Remove studies you don't need
studies_to_remove <- c("Vanheer 2025", "Lin  2017", "Myint 2014", "Yao 2023")

cleaned_existing_data <- existing_data %>%
  filter(!study %in% studies_to_remove)

cat("\nRemoved unnecessary studies:", paste(studies_to_remove, collapse = ", "), "\n")
cat("Remaining studies:", nrow(cleaned_existing_data), "rows\n")

# Prepare new data in same format as existing
new_data_formatted <- new_prevalence_data %>%
  select(
    study,
    location,
    latitude,
    longitude,
    date_start,
    date_end,
    prev
  ) %>%
  mutate(
    mean_age_yrs = NA,
    age_group = NA,
    year_s = paste(year(date_start), "to", year(date_end)),
    mid_date = date_start + (date_end - date_start)/2
  ) %>%
  # Reorder columns to match existing data
  select(names(cleaned_existing_data))

# Combine with existing data
complete_study_locations <- bind_rows(
  cleaned_existing_data,
  new_data_formatted
)

cat("\n=== FINAL DATASET ===\n")
cat("Total studies:", length(unique(complete_study_locations$study)), "\n")
cat("Total rows:", nrow(complete_study_locations), "\n")

# Show studies included
cat("\nStudies in final dataset:\n")
print(sort(unique(complete_study_locations$study)))

# Save the updated file
write_csv(complete_study_locations, "analysis/data-derived/study_locations_processed_complete.csv")

cat("\nSaved complete dataset to: analysis/data-derived/study_locations_processed_complete.csv\n")

# Now you can use this for your analysis
cat("\nNext step: Use 'complete_study_locations' for your enhanced regression analysis\n")
