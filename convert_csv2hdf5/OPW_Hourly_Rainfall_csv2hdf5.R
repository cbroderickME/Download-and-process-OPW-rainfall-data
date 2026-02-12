# ------------------------------------------------------------------------------
# This script processes hourly rainfall data from multiple stations spanning 2020 to 2025.
# It reads in station metadata and rainfall data, cleans and organizes the data,
# calculates statistics per station, and saves results both as individual CSV files 
# and within a consolidated HDF5 file.
#
# --- Summary of actions ---
# • Creates timestamp for output file naming.
# • Loads required R packages for data manipulation and HDF5 handling.
# • Creates necessary output directories if they do not exist.
# • Reads station metadata (station numbers and coordinates).
# • Creates or overwrites an HDF5 file to store processed data.
# • For each station:
#    - Filters rainfall data.
#    - Skips stations with only missing data (-99).
#    - Removes duplicate timestamps and fills missing hours with -99.
#    - Adds station metadata columns (Latitude, Longitude, station number).
#    - Calculates summary statistics:
#       - Start and end dates of the record.
#       - Maximum rainfall and the date it occurred.
#       - Largest gap (missing data streak) in the record.
#    - Saves individual station data as CSV.
#    - Writes station data into the HDF5 file under its own group.
# • After processing all stations:
#    - Aggregates and saves a summary CSV of statistics for all stations.
#    - Writes processed station metadata to the HDF5 file.
# • Demonstrates how to read back data from the HDF5 file (example).
# • Closes HDF5 resources at the end.
# ------------------------------------------------------------------------------

# ------------------------------ Step 1: Setup and Libraries ------------------------------
# Clear environment and set timezone
rm(list = ls())  
Sys.setenv(TZ = "GMT")

# Set current working directory and load necessary libraries
setwd("K:/Fluvial/Meteorological_Data/Rainfall/OPW")
wkdir <- getwd()  

timestamp <- format(Sys.time(), "%Y%m%d")  # Create timestamp for file names

max_time = (format(Sys.time(), "%Y"))

# Install and load packages (if not already installed)
library("dplyr")
library("lubridate") 
library("rhdf5") 
library("tibble") 
library("purrr") 
library("tidyr")

# Define the target directory
target_dir <- file.path(wkdir, "/csv")

# Delete everything in the directory if it exists
if (dir.exists(target_dir)) {
  unlink(target_dir, recursive = TRUE, force = TRUE)
}

# Recreate the directory (with parents if needed)
dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)

# Initialize results list
results <- list()

# ---------------------------- Step 2: Read in metadata table from climate database ----------------------------
# Load station metadata (station number, latitude, longitude)
stno_lst <- read.csv("K:/Fluvial/Meteorological_Data/Rainfall/OPW/HydroData_info/OPW_Pstats_metadata.csv", stringsAsFactors = FALSE) %>%
  as_tibble() %>%
  add_column("Processed" = NA)%>%
  rename(stno=ID)

# ---------------------------- Step 4: Create HDF5 File ----------------------------
fle_nme <- paste0(wkdir, "/hdf5/OPW_HrlyRain_2020_", max_time , "_Gen", timestamp, ".h5")
unlink(fle_nme)  # Remove if exists
h5createFile(fle_nme)


# ---------------------------- Step 5: Process Each Station ----------------------------
for (st_i in seq_along(stno_lst$stno)) {
  station_id <- stno_lst$stno[st_i]
  message("Processing station ", station_id, " (", st_i, " of ", nrow(stno_lst), ")")
  
  # Read station CSV
  file_path <- file.path(wkdir, "processed", paste0(station_id, "_hourly.csv"))
  if (!file.exists(file_path)) {
    warning("File not found: ", file_path, ". Skipping station.")
    next
  }
  
  Pid_con <- read.csv(file_path) %>%
    rename(rain = rainfall, date_time = time) %>%
    mutate(date_time = ymd_hms(date_time)) %>%
    distinct(date_time, .keep_all = TRUE)
  
  # Skip if all rainfall values are -99
  if (all(Pid_con$rain == -99, na.rm = TRUE)) {
    message("Skipping station ", station_id, " (all rainfall values are -99).")
    next
  }
  
  # Full hourly sequence
  full_time_seq <- tibble(date_time = seq(
    ymd_hms("2020-01-01 00:00:00"),
    Sys.time(),
    by = "hour"
  ))
  
  Pid_con <- full_time_seq %>%
    left_join(Pid_con, by = "date_time") %>%
    mutate(
      rain = ifelse(is.na(rain), -99, rain),
      Latitude = stno_lst$Latitude[st_i],
      Longitude = stno_lst$Longitude[st_i],
      stno = station_id,
      date_time = format(date_time, tz = "UTC", usetz = TRUE)
    )
  
  stno_lst$Processed[st_i] <- 'Y'
  
  # ---------------------- Statistics ----------------------
  start_date <- min(Pid_con$date_time)
  end_date <- max(Pid_con$date_time)
  max_rain <- max(Pid_con$rain, na.rm = TRUE)
  max_rain_date <- Pid_con %>% filter(rain == max_rain) %>% pull(date_time) %>% .[1]
  
  largest_gap <- Pid_con %>%
    mutate(missing = rain == -99) %>%
    mutate(group = cumsum(!missing)) %>%
    group_by(group) %>%
    summarize(
      start_time = first(date_time[missing]),
      end_time = last(date_time[missing]),
      duration = as.numeric(difftime(end_time, start_time, units = "hours")),
      .groups = 'drop'
    ) %>%
    filter(duration > 0) %>%
    slice_max(duration)
  
  if (nrow(largest_gap) == 0) {
    largest_gap <- tibble(start_time = NA, end_time = NA, duration = NA)
  }
  
  results[[length(results) + 1]] <- tibble(
    Station = station_id,
    Start_Date = start_date,
    End_Date = end_date,
    Max_Rain = max_rain,
    Max_Rain_Date = max_rain_date,
    Max_Missing_Duration = largest_gap$duration,
    Largest_Gap_Start = largest_gap$start_time,
    Largest_Gap_End = largest_gap$end_time
  )
  
  # ---------------------- Save Individual CSV ----------------------
  out_csv_dir <- file.path(target_dir, "Hourly_Rainfall_Individual_Stations")
  dir.create(out_csv_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(Pid_con, file = file.path(out_csv_dir,
                                      paste0("Hourly_Rainfall_", station_id, "_2020_", max_time, "_", timestamp, ".csv")),
            row.names = FALSE)
  
  # ---------------------- Write to HDF5 (column-wise) ----------------------
  h5createGroup(fle_nme, paste0(stno_lst$stno[st_i]))
  Pid_con$date_time <- as.character(Pid_con$date_time)
  h5write(Pid_con[, c("date_time", "rain", "Latitude", "Longitude")], file = fle_nme,
          paste0(stno_lst$stno[st_i], "/Hr_Rain"))
  
  
  message("Processed station: ", station_id)
}

# ---------------------------- Step 6: Finalize and Save Results ----------------------------
results_df <- bind_rows(results) %>% distinct()
results_file_path <- paste0(wkdir, "/csv/Hourly_Rainfall_Individual_Stations/Station_Summary_2020_", max_time ,"_", timestamp, ".csv")
write.csv(results_df, file = results_file_path, row.names = FALSE)
message("Results DataFrame saved to ", results_file_path)

# Save processed station metadata to the HDF5 file
meta <- stno_lst %>% filter(Processed == "Y") %>% select(stno, Latitude, Longitude)
h5write(meta, file = fle_nme, "Metadata")

message("Processing complete. HDF5 and CSV files saved.")
H5close()

# ---------------------------- Step 7: Read Data from HDF5 (Example) ----------------------------
# Example: Read data from HDF5 for the last processed station
station_data <- h5read(fle_nme, paste0(81632, "/Hr_Rain"))
station_df <- as.data.frame(station_data)
head(station_df)


