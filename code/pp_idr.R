## Postprocessing via IDR ##

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
# Get locations
loc_vec <- c("bon", "dra", "fpk", "gwn", 
             "psu", "sxf", "tbl")

# Training and test period
yr_tr <- 2017:2019
yr_te <- 2020

# Prediction variables
pred_vars_ls <- list("idr-ref" = "rest2",
                     "idr-ghi" = "ssrd")

#### IDR ####
# For-Loop over locations
for(temp_loc in loc_vec){
  #### Data ####
  # Load
  load(paste0(data_path, "data_", temp_loc, ".RData"))
  
  # train-test split
  data_tr <- data %>%
    filter(year(Time) %in% yr_tr)
  data_te <- data %>%
    filter(year(Time) %in% yr_te)
  
  #### For-Loop over predictor variables ####
  for(temp_meth in names(pred_vars_ls)){
    #### IDR ####
    # Console
    print(paste0(" --- Method: ", temp_meth, 
                 " --- Location: ", temp_loc, " ---"))
    
    # Fit IDR
    pred <- idr_pp(train = data_tr,
                   X = data_te,
                   pred_vars = pred_vars_ls[[temp_meth]],
                   idr_ls = list(n_sbg = 0))
    
    #### Save ####
    # Save
    save(file = paste0(data_path, temp_meth, "_", temp_loc, ".RData"),
         list = c("pred", "data_te"))
  }
}