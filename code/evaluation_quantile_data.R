## Generate data for quantile reliability diagrams ##

#### Housekeeping ####
rm(list=ls())
gc()

#### Settings ####
# Packages
library(dplyr)
library(lubridate)
library(tidyverse)

# Path of data
data_path <- "/data/"

# Path of R-functions
code_path <- "/code/"

# Load functions
setwd(code_path)
source(file = paste0(code_path, "functions_pp.R"))
source(file = paste0(code_path, "functions_quantile_rd.R"))

#### Initiation ####
# Prediction of bias with networks?
bqn_bias <- TRUE
drn_bias <- TRUE

# Get locations
loc_vec <- c("bon", "dra", "fpk", "gwn", 
             "psu", "sxf", "tbl")

# Postprocessing methods (local and global)
pp_local <- c("AnEn", "idr-ref", "idr-ghi")
pp_global <- c("drn", "bqn")
pp_vec <- c(pp_local, pp_global)

# Quantile levels
q_levels <- sort(unique(c((1:19)*5/100, 0.01, 0.02, 0.98, 0.99, c(1:11)/12)))

#### Create data frame ####
# Name of columns
col_vec <- c("date", "location", "quantile", "value", "model", "truth")

# Make data frame
df <- as.data.frame(matrix(nrow = 0,
                           ncol = length(col_vec)))
colnames(df) <- col_vec

# For-Loop over methods
for(temp_meth in pp_vec){
  # Check whether local or global
  if(is.element(temp_meth, pp_local)){
    # For-Loop over locations
    for(temp_loc in loc_vec){
      # Load data
      load(file = paste0(data_path, temp_meth, "_", temp_loc, ".RData"))
      
      # Create one data frame for each quantile
      fn_apply <- function(qu){
        # qu...Quantile level
        
        # Make data frame
        res <- as.data.frame(matrix(nrow = nrow(data_te),
                                    ncol = length(col_vec)))
        colnames(res) <- col_vec
        
        # Read out date, location, model and truth
        res[,"date"] <- data_te[["Time"]]
        res[,"location"] <- toupper(temp_loc)
        res[,"quantile"] <- qu
        res[,"model"] <- temp_meth
        res[,"truth"] <- data_te[["ghi"]]
        
        # Calculate quantiles
        if(temp_meth == "AnEn"){
          # Check if level corresponds to an ensemble member
          if(is.element(qu, (1:11)/12)){ res[,"value"] <- pred[["f"]][,12*qu] }
          else{ res[,"value"] <- apply(pred[["f"]], 1, function(x){ 
              quantile(x = x,
                       probs = qu,
                       type = 8) }) }
        }
        else if(grepl("idr", temp_meth, fixed = TRUE)){
          res[,"value"] <- as.vector(isodistrreg::qpred(predictions = pred[["pred_idr"]],
                                                        quantiles = qu))
        }
        
        # Output
        return(res)
      }
      
      # Apply on quantile levels
      df <- rbind(df, bind_rows(lapply(q_levels, fn_apply)))
    }
  }
  else if(is.element(temp_meth, pp_global)){
    # Get name of file
    if(((temp_meth == "drn") & drn_bias) |
       ((temp_meth == "bqn") & bqn_bias)){ temp_name <- paste0(temp_meth, "_bias") 
    }else{ temp_name <- temp_meth }
    
    # Load data
    load(file = paste0(data_path, temp_name, ".RData"))
    
    # Create one data frame for each quantile
    fn_apply <- function(qu){
      # qu...Quantile level
      
      # Make data frame
      res <- as.data.frame(matrix(nrow = nrow(data_te),
                                  ncol = length(col_vec)))
      colnames(res) <- col_vec
      
      # Read out date, location, model and truth
      res[,"date"] <- data_te[["Time"]]
      res[,"location"] <- toupper(data_te[["location"]])
      res[,"quantile"] <- qu
      res[,"model"] <- temp_meth
      res[,"truth"] <- data_te[["ghi"]]
      
      # Calculate quantiles
      if(temp_meth == "drn"){
        res[,"value"] <- crch::qtnorm(p = qu, 
                                      mean = pred[["f"]][,1], 
                                      sd = pred[["f"]][,2],
                                      left = 0)
      }
      else if(temp_meth == "bqn"){
        # Check whether quantile was already calculated
        if(is.element(qu, (1:99)/100)){
          # Read out quantile
          res[,"value"] <- pred[["f"]][,100*qu]
        }else{
          # Calculate quantiles based on coefficients (Use 1 - qu due to change in orientation)
          res[,"value"] <- as.vector(bern_quants(alpha = pred[["alpha"]],
                                                 q_levels = 1 - qu))
          
          # Transform bias to target variable
          if(bqn_bias){ res[,"value"] <- pmax(data_te[,"ssrd"] - res[,"value"], 0) }
        }
      }
      
      # Output
      return(res)
    }
    
    # Apply on quantile levels
    df <- rbind(df, bind_rows(lapply(q_levels, fn_apply)))
  }
  rm(pred, data_te)
}

# Save data frame
save(file = paste0(data_path, "data_quantile_rd.RData"),
     list = c("df"))
