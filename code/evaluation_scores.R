## summarize evaluation measures ##

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

# Quantile levels for Brier score and pinball loss
q_levels <- c(11/12, 0.9, 0.95, 0.99)
q_names <- round(100*q_levels, 2)

# Threshold levels for Brier score
t_levels <- c(250, 500, 700, 750, 800)

# Significance niveau of prediction intervals
pi_alpha <- 2/12
pi_levels <- c(pi_alpha/2, 1 - pi_alpha/2)

# Evaluation measures
sr_vec <- c("crps", "logs", "me", "mae", "rmse", "is", "lgt", "cov",
            paste0("bs_p", q_names), paste0("bs_t", t_levels), paste0("pinball", q_names))

# Column names of locationwise data frame
col_vec <- c("method", "location", "n", "runtime_est", "runtime_pred", sr_vec)

# Make data frame
df_scores <- as.data.frame(matrix(ncol = length(col_vec),
                                  nrow = length(pp_vec)*length(loc_vec)))
colnames(df_scores) <- col_vec

# List for PIT values and quantiles
pit_ls <- t_q_ls <- list()

#### Calculate locationwise quantiles ####
# For-Loop over locations
for(temp_loc in loc_vec){
  # Load data
  load(file = paste0(data_path, "data_", temp_loc, ".RData"))
  
  # Calculate thresholds
  t_q_ls[[temp_loc]] <- quantile(x = data[["ghi"]],
                                 probs = q_levels,
                                 type = 8)
}

#### Calculate measures ####
# For-Loop over methods
for(temp_meth in pp_vec){
  # Check whether local or global
  if(is.element(temp_meth, pp_local)){
    # For-Loop over locations
    for(temp_loc in loc_vec){
      # Load data
      load(file = paste0(data_path, temp_meth, "_", temp_loc, ".RData"))
      
      # Get index
      i_df <- which(temp_loc == loc_vec) + (which(temp_meth == pp_vec) - 1)*length(loc_vec)
      
      # Write in data frame
      df_scores[i_df, "method"] <- temp_meth
      df_scores[i_df, "location"] <- temp_loc
      df_scores[i_df, "n"] <- nrow(data_te)
      df_scores[i_df, "runtime_est"] <- pred[["runtime_est"]]
      df_scores[i_df, "runtime_pred"] <- pred[["runtime_pred"]]
      df_scores[i_df, "crps"] <- mean(pred[["scores_pp"]][["crps"]])
      df_scores[i_df, "me"] <- mean(pred[["scores_pp"]][["e_md"]])
      df_scores[i_df, "mae"] <- mean(abs(pred[["scores_pp"]][["e_md"]]))
      df_scores[i_df, "rmse"] <- sqrt(mean(pred[["scores_pp"]][["e_me"]]^2))
      
      # Prediction interval refers to range of a 11-member ensemble (AnEn)
      df_scores[i_df, "cov"] <- fn_cover(x = pred[["scores_pp"]][["pit"]],
                                         alpha = pi_alpha)

      # Get lower and upper bounds of prediction interval
      if(grepl("idr", temp_meth, fixed = TRUE)){
        # Calculate quantiles
        pi_bounds <- isodistrreg::qpred(predictions = pred[["pred_idr"]], 
                                      quantiles = pi_levels)
      }
      else if(is.element(temp_meth, c("AnEn"))){
        # Check levels
        if(pi_alpha == 2/12){ pi_bounds <- pred[["f"]][,c(1, 11)]
        }else{ pi_bounds0 <- t(apply(pred[["f"]], 1, function(x){
          quantile(x = x,
                   probs = pi_levels,
                   type = 8)})) }
      }
      
      # Calculate prediction interval length
      df_scores[i_df, "lgt"] <- mean(pi_bounds[,2] - pi_bounds[,1])
      
      # Calculate interval score
      df_scores[i_df, "is"] <- mean(interval_score(l = pi_bounds[,1],
                                                   u = pi_bounds[,2],
                                                   y = data_te[["ghi"]]))
      
      # Different Brier score calculations
      if(grepl("idr", temp_meth, fixed = TRUE)){
        # Function for Brier score depending on quantile level x
        fn_bs <- function(x){
          mean(isodistrreg::bscore(predictions = pred[["pred_idr"]], 
                                   threshold = x, 
                                   y = data_te[["ghi"]]))
        }
      }
      else if(is.element(temp_meth, c("AnEn"))){
        fn_bs <- function(x){ mean(brier_score(f = pred[["f"]],
                                               y = data_te[["ghi"]],
                                               t = x,
                                               distr = "ens")) }
      }
      
      # Different pinball loss calculations
      if(grepl("idr", temp_meth, fixed = TRUE)){
        # Function for pinball loss (in package scaled by factor 2!)
        fn_pinball <- function(x){
          mean(isodistrreg::qscore(predictions = pred[["pred_idr"]], 
                                   quantiles = x, 
                                   y = data_te[["ghi"]]))/2
        }
      }
      else if(is.element(temp_meth, c("AnEn"))){
        fn_pinball <- function(x){ 
          # Check levels
          if(x == 11/12){ res <- pinball_loss(f = pred[["f"]][,11],
                                              y = data_te[["ghi"]],
                                              alpha = x,
                                              distr = "q") }
          else{ res <- pinball_loss(f = pred[["f"]],
                                    y = data_te[["ghi"]],
                                    alpha = x,
                                    distr = "ens") }
          
          return(mean(res)) 
        }
      }
      
      # Calculate Brier scores and pinball loss
      df_scores[i_df, paste0("bs_p", q_names)] <- sapply(t_q_ls[[temp_loc]], fn_bs)
      df_scores[i_df, paste0("bs_t", t_levels)] <- sapply(t_levels, fn_bs)
      df_scores[i_df, paste0("pinball", q_names)] <- sapply(q_levels, fn_pinball)
      
      # Read out PIT values
      pit_ls[[paste0(temp_meth, "_", temp_loc)]] <- pred[["scores_pp"]][["pit"]]
    }
  }
  else if(is.element(temp_meth, pp_global)){
    # Get name of file
    if(((temp_meth == "drn") & drn_bias) |
       ((temp_meth == "bqn") & bqn_bias)){ temp_name <- paste0(temp_meth, "_bias") 
    }else{ temp_name <- temp_meth }
    
    # Load data
    load(file = paste0(data_path, temp_name, ".RData"))
    
    # For-Loop over locations
    for(temp_loc in loc_vec){
      # Get indices of location
      i_loc <- which(data_te[["location"]] == temp_loc)
      
      # Subset of scores
      scores_loc <- pred[["scores_pp"]][i_loc,]

      # Get index
      i_df <- which(temp_loc == loc_vec) + (which(temp_meth == pp_vec) - 1)*length(loc_vec)
      
      # Write in data frame
      df_scores[i_df, "method"] <- temp_meth
      df_scores[i_df, "location"] <- temp_loc
      df_scores[i_df, "n"] <- length(i_loc)
      df_scores[i_df, "runtime_est"] <- pred[["runtime_est"]]
      df_scores[i_df, "runtime_pred"] <- pred[["runtime_pred"]]*df_scores[i_df, "n"]/nrow(data_te)
      df_scores[i_df, "crps"] <- mean(scores_loc[["crps"]])
      df_scores[i_df, "me"] <- mean(scores_loc[["e_md"]])
      df_scores[i_df, "mae"] <- mean(abs(scores_loc[["e_md"]]))
      df_scores[i_df, "rmse"] <- sqrt(mean(scores_loc[["e_me"]]^2))
      
      # Prediction interval refers to range of a 11-member ensemble (AnEn)
      df_scores[i_df, "cov"] <- fn_cover(x = scores_loc[["pit"]],
                                         alpha = pi_alpha)
      
      # Locationwise forecasts and observations
      f_loc <- pred[["f"]][i_loc,]
      y_loc <- data_te[i_loc, "ghi"]
      
      # Get lower and upper bounds of prediction interval
      if(temp_meth == "drn"){
        # Get correct distribution
        if(drn_bias){ temp_distr <- paste0("t", pred[["nn_ls"]][["distr"]]) 
        }else{ temp_distr <- pred[["nn_ls"]][["distr"]] }
        
        # Calculate quantiles
        if(temp_distr == "tnorm"){ pi_bounds <- 
            sapply(pi_levels, function(x) crch::qtnorm(p = x,
                                                       mean = f_loc[,1],
                                                       sd = f_loc[,2],
                                                       left = 0)) 
        }else if(temp_distr == "tlogis"){ pi_bounds <- 
          sapply(pi_levels, function(x) crch::qtlogis(p = x,
                                                      location = f_loc[,1],
                                                      scale = f_loc[,2],
                                                      left = 0)) }
      }
      else if(temp_meth == "bqn"){
        # Sum up calculated quantiles of bias (Sum of basis at quantiles times coefficients)
        pi_bounds <- bern_quants(alpha = pred[["alpha"]][i_loc,],
                                 q_levels = pi_levels)
        
        # Add to forecasts for bias
        if(bqn_bias){
          # Add up forecast and bias while ensuring positivity
          pi_bounds <- pmax(data_te[["ssrd"]][i_loc] - pi_bounds, 0) 
          
          # Orientation has changed
          pi_bounds <- t(apply(pi_bounds, 1, rev))
        }
      }
      
      # Calculate prediction interval length
      df_scores[i_df, "lgt"] <- mean(pi_bounds[,2] - pi_bounds[,1])
      
      # Calculate interval score
      df_scores[i_df, "is"] <- mean(interval_score(l = pi_bounds[,1],
                                                   u = pi_bounds[,2],
                                                   y = y_loc))
      
      # Different Brier score calculations
      if(temp_meth == "drn"){
        # Get correct distribution
        if(drn_bias){ temp_distr <- paste0("t", pred[["nn_ls"]][["distr"]]) 
        }else{ temp_distr <- pred[["nn_ls"]][["distr"]] }
        
        fn_bs <- function(x){
          mean(brier_score(f = f_loc,
                           y = y_loc,
                           t = x,
                           distr = temp_distr,
                           t_distr = 0)) }
      }
      else if(temp_meth == "bqn"){
        fn_bs <- function(x){ mean(brier_score(f = f_loc,
                                               y = y_loc,
                                               t = x,
                                               distr = "ens")) }
      }
      
      # Different pinball loss calculations
      if(temp_meth == "drn"){
        # Get correct distribution
        if(drn_bias){ temp_distr <- paste0("t", pred[["nn_ls"]][["distr"]]) 
        }else{ temp_distr <- pred[["nn_ls"]][["distr"]] }
        
        fn_pinball <- function(x){ mean(pinball_loss(f = f_loc,
                                                     y = y_loc,
                                                     alpha = x,
                                                     distr = temp_distr,
                                                     t_distr = 0)) }
      }
      else if(temp_meth == "bqn"){
        fn_pinball <- function(x){ 
          # Read out quantile forecasts (based on assumption of 99 equidistant quantiles)
          q_loc <- f_loc[,100*x]

          # Calculate scores
          mean(pinball_loss(f = q_loc,
                            y = y_loc,
                            alpha = x,
                            distr = "q"))
        }
      }
      
      # Calculate Brier scores and pinball loss
      df_scores[i_df, paste0("bs_p", q_names)] <- sapply(t_q_ls[[temp_loc]], fn_bs)
      df_scores[i_df, paste0("bs_t", t_levels)] <- sapply(t_levels, fn_bs)
      df_scores[i_df, paste0("pinball", q_names)] <- sapply(q_levels, fn_pinball)
      
      # Read out PIT values
      pit_ls[[paste0(temp_meth, "_", temp_loc)]] <- scores_loc[["pit"]]
    }
  }
}

#### Save ####
# Save
save(file = paste0(data_path, "scores.RData"),
     list = c("df_scores", "pit_ls"))