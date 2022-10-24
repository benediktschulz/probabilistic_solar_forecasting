## Postprocessing via AnEn ##

# Here, we apply the postprocessing method applied in Yang et al. (2022).
# The code is taken from the file "AnEn.R" in the corresponding supplemental material.

#### Housekeeping ####
rm(list=ls())
gc()

#### Settings ####
# Packages
library(dplyr)
library(lubridate)
library(RANN)

# Path of data
data_path <- "/data/"

# Path of R-functions
code_path <- "/code/"

# Path of original data (Yang et al., 2022)
paper_path <- paste0(data_path, "original_data/")

# Load functions
setwd(code_path)
source(file = paste0(code_path, "functions_eval.R"))

#### Code from Paper ####
#################################################################################
# Inputs
#################################################################################
# dir.data <- "/Volumes/Macintosh Research/Data Article" # highest-level directory 
# dir.processed <- file.path(dir.data, "Processed") # directory for ECMWF data
# source("/Users/dyang/Dropbox/Working papers/ECMWF3/Code/functions.R") # source functions
yr.tr <- c(2017, 2018, 2019) # training years
yr.te <- c(2020) # verification years
stn <- c("bon", "dra", "fpk", "gwn", "psu", "sxf", "tbl")
Q <- c(0.025, seq(0.1, 0.90, by = 0.1), 0.975) # a set of quantiles to evaluate
m <- length(Q) # number of ensemble members
#################################################################################

# read metadata of SURFRAD stations
load(paste0(paper_path, "SURFRAD.loc.RData"))
loc <- SURFRAD.loc
# reorder rows of "loc", to follow the order of the stations
loc <- loc[match(stn, loc$stn),]

# loop over the stations
# set a few empty lists to hold the forecasts and observations
f <- list()
pb = txtProgressBar(min = 0, max = nrow(loc), initial = 0, style = 3) 
for(stn.ind in 1:nrow(loc))
{
  # read previously arranged obs and fcst datasets
  setwd(paper_path)
  load(paste0(stn[stn.ind], "_obs.RData")) # load observations
  load(paste0(stn[stn.ind], "_fcst.RData")) # load forecasts
  
  # combine the two tibbles in to one
  # change from UTC to local time
  # remove the first and last day, which has incomplete data due to time zone change
  # compute forecast clear-sky index (if index > 1.2, cap it at 1.2)
  # convert temperatures to degrees
  # calculate relative humidity from dew point and ambient temperature
  # compute bias in the forecasts
  data <- obs %>%
    left_join(., fcst, by = "Time") %>%
    mutate(Time = Time + loc$tz[stn.ind]*3600) %>%
    filter(date(Time) >= date("2017-01-02") & date(Time) <= date("2020-12-30")) %>%
    mutate(kappa = ifelse(ssrd/rest2 > 1.2, 1.2, ssrd/rest2)) %>%
    mutate(t2m = t2m - 273.15) %>%
    mutate(d2m = d2m - 273.15) %>%
    mutate(rh = 10^(7.591386*(d2m/(d2m+240.7263)-t2m/(t2m+240.7263))))
  
  # define feature columns
  feature_cols <- c("kappa", "zen", "t2m", "sp", "rh")
  
  # train-test split
  # filter off zen > 85 points
  data.tr <- data %>%
    filter(year(Time) %in% yr.tr) %>%
    filter(zen < 85)
  data_te <- data %>%
    filter(year(Time) %in% yr.te) %>%
    filter(zen < 85)
  Xtr <- data.matrix(data.tr[,feature_cols])
  Xte  <- data.matrix(data_te[,feature_cols])
  
  # nomalization of training data, and re-scale the test data with the same values
  Xtr <- scale(Xtr)
  Xtr[is.nan(Xtr)] <- 0 # python uses 0 by default, in R, this is done manually
  Xte = scale(Xte, center=attr(Xtr, "scaled:center"), scale=attr(Xtr, "scaled:scale"))
  Xte[is.nan(Xte)] <- 0
  
  # Take time
  start_tm <- Sys.time()
  
  # kd-tree algorithm, choose 21 analogs
  analogs_kdtree <- RANN::nn2(Xtr, query = Xte, k = m, treetype = "kd", searchtype = "standard")
  analogs_kdtree <- analogs_kdtree$nn.idx
  
  # select the corresponding observations
  Ytr <- data.tr$ghi/data.tr$rest2 # observed historical clear-sky index
  kappa_pred <- apply(analogs_kdtree, 2, function(x) Ytr[x])
  ghi_pred <- apply(kappa_pred, 2, function(x) x*data_te$rest2)
  f <- t(apply(ghi_pred, 1, sort)) # sort the quantiles
  
  # Take time
  end_tm <- Sys.time()
  
  # Time needed
  runtime_est <- as.numeric(difftime(end_tm, start_tm, units = "mins"))
  
  # Evaluation
  scores_pp <- fn_scores_ens(ens = f,
                             y = data_te$ghi)
  
  # Transform ranks to uPIT-values
  scores_pp[["pit"]] <- scores_pp[["rank"]]/(ncol(f) + 1) - 
    runif(n = nrow(f),
          min = 0,
          max = 1/(ncol(f) + 1))
  scores_pp[["rank"]] <- NULL
  
  # Save data in list
  pred <- list(f = f,
               scores_pp = scores_pp,
               runtime_est = runtime_est,
               runtime_pred = NA)
  
  # Save analog forecasts and data
  save(file = paste0(data_path, "AnEn_", stn[stn.ind], ".RData"),
       list = c("pred", "data_te"))
  
  setTxtProgressBar(pb,stn.ind)
}
close(pb)