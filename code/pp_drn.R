## Postprocessing via DRN ##

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

# Load functions
setwd(code_path)
source(file = paste0(code_path, "functions_pp.R"))

#### Initiation ####
# Predict bias?
log_bias <- TRUE

# Get locations
loc_vec <- c("bon", "dra", "fpk", "gwn", 
             "psu", "sxf", "tbl")

# Training and test period
yr_tr <- 2017:2019
yr_te <- 2020

# Prediction variables
pred_vars <- c("kappa", "zen", "t2m", "sp", "rh", "location")

#### Data ####
# Load
load(paste0(data_path, "data_total.RData"))

# train-test split
data_tr <- data %>%
  filter(year(Time) %in% yr_tr)
data_te <- data %>%
  filter(year(Time) %in% yr_te)
rm(data)

#### DRN ####
# Function for prediction
if(log_bias){ fn_drn <- drn_pp_bias 
}else{ fn_drn <- drn_pp }

# Call DRN
pred <- fn_drn(train = data_tr,
               X = data_te,
               i_valid = which(year(data_tr$Time) == 2019),
               loc_id_vec = loc_vec,
               pred_vars = pred_vars)

#### Save ####
# File name
temp_name <- "drn"

# Bias prediction?
if(log_bias){ temp_name <- paste0(temp_name, "_bias")}

# Save
save(file = paste0(data_path, temp_name, ".RData"),
     list = c("pred", "data_te"))