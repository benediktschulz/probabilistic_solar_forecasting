## Data preprocessing ##

# Here, we apply the preprocessing steps that were performed by Yang et al. (2022).
# The code is taken from the file "AnEn.R" in the corresponding supplemental material.
# The underlying data is also taken from the supplemental material.

#### Housekeeping ####
rm(list=ls())
gc()

#### Settings ####
# Packages
library(dplyr)
library(lubridate)

# Path of data
data_path <- "/data/"

# Path of R-functions
code_path <- "/code/"

# Path of original data (Yang et al., 2022)
paper_path <- paste0(data_path, "original_data/")

# Load functions
setwd(code_path)

#### Initiation ####
# Get list of files
file_names <- list.files(paper_path)

# Get locations
loc_vec <- unique(substr(file_names, 1, 3))

# Difference d to UTC (UTC + d = local)
d_utc <- c("bon" = 6,
           "dra" = 8,
           "fpk" = 7,
           "gwn" = 6,
           "psu" = 5,
           "sxf" = 6,
           "tbl" = 7)

#### Data processing ####
# List for total data
data_ls <- list()

# For-Loop over locations
for(temp_loc in loc_vec){
  # Load forecast and observational data
  load(paste0(paper_path, temp_loc, "_fcst.RData"))
  load(paste0(paper_path, temp_loc, "_obs.RData"))

  # Combine two data sets
  data <- obs %>%
    left_join(., fcst, by = "Time") %>% # Combine two sets
    mutate(Time = Time - hours(d_utc[temp_loc])) %>% # Local time
    filter(date(Time) >= date("2017-01-02") & date(Time) <= date("2020-12-30")) %>% # Skip first and last dates
    mutate(kappa = ifelse(ssrd/rest2 > 1.2, 1.2, ssrd/rest2)) %>% # Forecast clear-sky index (if index > 1.2, cap it at 1.2)
    mutate(t2m = t2m - 273.15) %>% # Convert to degree
    mutate(d2m = d2m - 273.15) %>% # Convert to degree
    mutate(rh = 10^(7.591386*(d2m/(d2m+240.7263)-t2m/(t2m+240.7263)))) %>% # Relative humidity
    mutate(bias = ssrd - ghi)  # Forecast bias

  # Filter off zen >= 85 points
  data <- data %>% filter(zen < 85)

  # Save in list
  data_ls[[temp_loc]] <- cbind("location" = temp_loc, data)
  
  # Save data
  save(file = paste0(data_path, "data_", temp_loc, ".RData"),
       list = c("data"))
}

# Bind data
data <- bind_rows(data_ls)

# Sort w.r.t. date
data <- data[order(data$Time),]

# Reset row names
row.names(data) <- 1:nrow(data)

# Save data
save(file = paste0(data_path, "data_total.RData"),
     list = c("data"))
