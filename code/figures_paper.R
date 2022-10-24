## Figures Paper ##

#### Housekeeping ####
rm(list=ls())
gc()

#### Settings ####
# Packages
library(dplyr)
library(lubridate)
library(ggplot2)
library(xtable)
library(tidyverse)
library(reliabilitydiag)

# Path of data
data_path <- "/data/"

# Path of R-functions
code_path <- "/code/"

# Path of Figures
fig_path <- "/figures/"

# Load functions
setwd(code_path)
source(file = paste0(code_path, "functions_pp.R"))
source(file = paste0(code_path, "replicate_dw_reliability_functions.R"))

#### Initiation ####
# Prediction of bias with networks?
bqn_bias <- TRUE
drn_bias <- TRUE

# Get locations
loc_vec <- c("bon", "dra", "fpk", "gwn", 
             "psu", "sxf", "tbl")

# Location to focus on
loc_plot <- "bon"

# Postprocessing methods (local and global)
pp_local <- c("idr-ref", "AnEn")
pp_global <- c("drn", "bqn")
pp_vec <- c(pp_local, pp_global)

# Order of methods for table
pp_order <- pp_vec

# Name of postprocessing methods
pp_labels <- sapply(pp_order, function(x) 
  if(x == "AnEn"){ return(x) }
  else if(x == "idr-ref"){ return("CSD-IDR") }
  else{ return(toupper(x)) })

#### Load data ####
# Save
load(file = paste0(data_path, "scores.RData"))

#### Forecast illustration ####
# Initialization time
init_tm <- str2time(c("2020071000"))

# Quantile levels
q_vec <- (1:11)/12

# Quantiles to be highlighted in blue
highlight <- (1:11)/12

# Number of quantile levels
n_q <- length(q_vec)

# Columns of data frame
col_vec <- c("yday", "hour", "quantile", "value", "model", "truth")

# Day of the year to plot
yday_plot <- yday(init_tm)

# Make data frame
df <- data.frame(matrix(nrow = 0,
                        ncol = length(col_vec)))
colnames(df) <- col_vec

# For-Loop over methods
for(temp_meth in pp_vec){
  # Differentiate local and global postprocessing
  if(is.element(temp_meth, pp_global)){
    # Get name of file
    if(((temp_meth == "drn") & drn_bias) |
       ((temp_meth == "bqn") & bqn_bias)){ temp_name <- paste0(temp_meth, "_bias") 
    }else{ temp_name <- temp_meth }
    
    # Load data
    load(file = paste0(data_path, temp_name, ".RData"))
    
    # Get indices of dates of interest
    i_sub <- which((is.element(yday(data_te[["Time"]]), yday_plot)) &
                     (data_te[["location"]] == loc_plot))
  }
  else{ 
    # Load data
    load(file = paste0(data_path, temp_meth, "_", loc_plot, ".RData")) 
    
    # Get indices of dates of interest
    i_sub <- which(is.element(yday(data_te[["Time"]]), yday_plot))
  }
  
  # Subset of data
  data_sub <- data_te[i_sub,]
  
  # Initiate forecast matrix
  f_sub <- matrix(nrow = nrow(data_sub),
                  ncol = length(q_vec))
  
  # Get quantile forecasts
  if(grepl("idr", temp_meth, fixed = TRUE)){
    # Check which quantiles have been calculated
    i_q <- is.element(q_vec, pred[["q_levels"]])
    
    # Read out given quantiles
    if(any(i_q)){ f_sub[,i_q] <- 
      sapply(q_vec[i_q], function(x) pred[["f"]][i_sub, which(pred[["q_levels"]] == x)]) }
    
    # Calculate remaining quantiles
    if(any(!i_q)){ f_sub[,!i_q] <- 
      isodistrreg::qpred(predictions = pred[["pred_idr"]],
                         quantiles = q_vec[!i_q])[i_sub,] }
  }
  else if(temp_meth == "bqn"){
    # Check which quantiles have been calculated
    i_q <- is.element(q_vec, pred[["q_levels"]])
    
    # Read out given quantiles
    if(any(i_q)){ f_sub[,i_q] <- 
      sapply(q_vec[i_q], function(x) pred[["f"]][i_sub, which(pred[["q_levels"]] == x)]) }
    
    # Calculate remaining quantiles
    if(any(!i_q)){ f_sub[,!i_q] <- pmax(data_sub[["ssrd"]] - 
                                          bern_quants(alpha = pred[["alpha"]][i_sub,],
                                                      q_levels = 1 - q_vec[!i_q]), 0) }
  }
  else if(temp_meth == "drn"){
    # DRN forecasts drawn from distribution
    f_sub <- sapply(q_vec, function(x) crch::qtnorm(p = x,
                                                    mean = pred[["f"]][i_sub, 1],
                                                    sd = pred[["f"]][i_sub, 2],
                                                    left = 0))
  }
  else if(temp_meth == "AnEn"){
    # Check which quantiles have been calculated
    i_q <- is.element(q_vec, (1:11)/12)
    
    # Read out given quantiles
    if(any(i_q)){ f_sub[,i_q] <- 
      sapply(q_vec[i_q], function(x) pred[["f"]][i_sub, which((1:11)/12 == x)]) }
    
    # Calculate remaining quantiles
    if(any(!i_q)){ f_sub[,!i_q] <- 
      t(sapply(i_sub, function(i) quantile(x = pred[["f"]][i,],
                                           probs = q_vec[!i_q],
                                           type = 8))) }
  }
  
  # Number of forecasts
  n_fc <- nrow(data_sub)
  
  # Indices of method
  i_df <- 1:(n_fc*n_q) + (n_fc*n_q)*(which(pp_vec == temp_meth) - 1)
  
  # Read out variables
  df[i_df, "model"] <- temp_meth
  df[i_df, "yday"] <- rep(x = yday(data_sub[["Time"]]),
                          each = n_q)
  df[i_df, "hour"] <- rep(x = hour(data_sub[["Time"]]),
                          each = n_q)
  df[i_df, "truth"] <- rep(x = data_sub[["ghi"]],
                           each = n_q)
  df[i_df, "quantile"] <- rep(x = q_vec,
                              times = n_fc)
  df[i_df, "value"] <- as.vector(t(f_sub))
}

# Order methods
df[["model"]] <- factor(df[["model"]],
                        levels = pp_order,
                        labels = pp_labels)

# Get dates
init_labels <- sapply(init_tm, function(t) 
  paste0(month.name[month(t)], " ", day(t), ", ", year(t)))

# Rename day of the year
df[["yday"]] <- factor(df[["yday"]],
                       levels = yday_plot,
                       labels = init_labels)

# Quantiles to be highlighted (and not)
df1 <- df %>% filter(quantile %in% highlight)
df2 <- df %>% filter(!quantile %in% highlight)

# Quantiles used for drawing boxes
df_box <- df %>%
  filter(quantile %in% c(min(q_vec), 0.5, max(q_vec))) %>% 
  pivot_wider(names_from = quantile, names_prefix = "value.", values_from = value)

# Rename columns of min and max
colnames(df_box)[colnames(df_box) == paste0("value.", min(q_vec))] <- "value.min"
colnames(df_box)[colnames(df_box) == paste0("value.", max(q_vec))] <- "value.max"

# Start plot
pdf_plot <- ggplot(df1, aes(x = hour)) 

# Panel depending on months and methods
if(length(yday_plot) == 1){
  pdf_plot <- pdf_plot +
    facet_wrap("model",
               ncol = 4,
               scales = "fixed")
}else{
  pdf_plot <- pdf_plot +
    facet_grid(rows = vars(model), 
               cols = vars(yday), 
               scales = "free_x")
}

# Continue plotting
pdf_plot <- 
  pdf_plot +
  geom_crossbar(data = df_box, aes(y = value.0.5, 
                                   ymin = value.min, 
                                   ymax = value.max), 
                fatten = 1, width = 0.4, size = 0.3, colour = "azure4", fill = "gray", alpha = 0.4) +
  geom_segment(data = df2, aes(x = hour - 0.2, xend = hour + 0.2, y = value, yend = value), 
               color = "azure4", size = 0.25) + 
  geom_segment(data = df1, aes(x = hour - 0.25, xend = hour + 0.25, y = value, yend = value), 
               color = "deepskyblue4", size = 0.4, lineend = "round") + 
  geom_line(aes(y = truth, col = 'darkred'), size = 0.5) +
  geom_segment(data = df1, aes(x = hour - 0.4, xend = hour + 0.4, y = truth, yend = truth), 
               color = "darkred", size = 0.4, lineend = "round") + 
  scale_color_identity(name = NULL,
                       breaks = c("darkred"),
                       labels = c("Truth"),
                       guide = "legend") +
  scale_x_continuous(breaks = c(6, 9, 12, 15, 18)) +
  xlab("Hour of the Day") +
  ylab(expression(paste("GHI [W/", m^2, "]"))) +
  theme_bw(base_size = 11) +
  theme(panel.grid.major = element_line(size = 0.05), 
        panel.grid.minor = element_blank(),
        legend.position = c(0, 1),
        legend.justification = c(0, 1),
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12))

# Save as PDF
ggsave(filename = paste0("fc_illustration_", loc_plot, ".pdf"),
       path = fig_path,
       plot = pdf_plot,
       width = 20,
       height = 9,
       scale = 0.6)

#### Evaluation metric tables ####
# Evaluation measures
sr_vec <- c("crps", "is", "lgt", "cov", "bs_p95", 
            "bs_t250", "bs_t300", "bs_t500", "bs_t750") # Pinball via Q-reldiag

# For-Loop over evaluation measures
for(temp_sr in sr_vec){
  # Create data frame
  df_table <- as.data.frame(matrix(nrow = length(pp_vec),
                                   ncol = length(loc_vec)))
  colnames(df_table) <- toupper(loc_vec)
  
  # Set row names
  row.names(df_table) <- pp_labels

  # For-Loop over methods
  for(temp_meth in pp_order){ for(temp_loc in loc_vec){
    # Index of score
    i_score <- which((df_scores[["method"]] == temp_meth) &
                       (df_scores[["location"]] == temp_loc))
    
    # Write in table
    df_table[which(pp_order == temp_meth), which(loc_vec == temp_loc)] <- 
      df_scores[i_score, temp_sr]
  }}
  
  # Digits
  if(grepl("bs", temp_sr, fixed = TRUE)){ temp_digits <- 3 
  }else if(is.element(temp_sr, c("crps", "lgt"))){ 
    temp_digits <- 1 
  }else{ temp_digits <- 2 }
  
  # Make xtable object
  df_table <- xtable(df_table, 
                     type = "latex", 
                     digits = temp_digits)
  
  # Columns
  align(df_table) <- paste0("l", paste0(rep("c", length(loc_vec)), collapse = ""))
  
  # Header
  addtorow <- list()
  addtorow$pos <- list(0)
  addtorow$command <- paste0("Method & ", paste0(toupper(loc_vec), collapse = " & "), " \\\\\n")
  
  # Horizontal lines to include
  hline.after <- c(-1, 0, nrow(df_table))
  
  # Make Latex table
  print(file = paste0(fig_path, "table_", temp_sr, ".txt"),
        x = df_table, 
        include.colnames = FALSE,
        add.to.row = addtorow,
        hline.after = hline.after,
        booktabs = TRUE,
        # only.contents = TRUE,
        format.args = list(big.mark = ",", 
                           decimal.mark = "."))
}

#### Quantile score decomposition tables ####
# Quantile levels
q_levels <- c(1, 11)/12

# Number of resamples (Random small number for which it works)
n_resamples <- 29 # Scores are independent of resamples

# Load data frame
load(file = paste0(data_path, "data_quantile_rd.RData"))

# For-Loop over qauntile levels
for(q_level in q_levels){
  # Subset of data
  df_sub <- subset(df, (quantile == q_level) & 
                     is.element(model, pp_vec))
  
  # Name of file for quantile reliability diagram data
  file_name <- paste0(data_path, "qreldiag_", n_resamples, "_", paste0(pp_vec, collapse = "_"), 
                      "_q", round(100*q_level, 2), ".RData")
  
  # Check if file exists
  if(file.exists(file_name)){ load(file = file_name)
  }else{
    # List for data
    ls_reldiag <- list()
  
    # For-Loop over locations
    for(temp_loc in loc_vec){
      # Subset of location
      df_loc <- subset(df_sub, location == toupper(temp_loc))
  
      # Data frame for plotting
      df_reldiag <- df_loc %>%
        group_by(model, quantile) %>%
        summarize(reldiag(value, truth, alpha = unique(quantile), 
                          n_resamples = n_resamples, digits = 2, ties = "primary"),
                  .groups = "keep") %>%
        mutate(across(c(x_rc, lower, upper), ~ pmax(., 0))) # set negative values to zero
  
      # Save in list
      ls_reldiag[[temp_loc]] <- cbind(df_reldiag, "location" = temp_loc)
    }
    rm(df_reldiag)
  
    # Save data
    save(file = file_name,
         list = c("ls_reldiag"))
  }
  
  # Create tables for pinball loss and decomposition components
  pinball_tbl <- mcb_tbl <- dsc_tbl <- 
    as.data.frame(matrix(nrow = length(pp_vec),
                         ncol = length(loc_vec)))
  
  # Set column and row names
  colnames(pinball_tbl) <- colnames(mcb_tbl) <- colnames(dsc_tbl) <- toupper(loc_vec)
  row.names(pinball_tbl) <- row.names(mcb_tbl) <- row.names(dsc_tbl) <- pp_labels
  
  # For-Loop over methods and locations
  for(temp_loc in toupper(loc_vec)){ for(temp_meth in pp_labels){
    # Get subset
    df_reldiag <- ls_reldiag[[tolower(temp_loc)]]
    
    # Get subset
    df_tbl <- subset(df_reldiag, (location == tolower(temp_loc)) &
                       (model == pp_order[which(pp_labels == temp_meth)]))
    
    # Get scores (which are unique)
    pinball_tbl[temp_meth, temp_loc] <- unique(df_tbl[["score"]])
    mcb_tbl[temp_meth, temp_loc] <- unique(df_tbl[["mcb"]])
    dsc_tbl[temp_meth, temp_loc] <- unique(df_tbl[["dsc"]])
  }}
  
  # For-Loop over three tables
  for(temp in c("pinball", "mcb", "dsc")){
    # Make xtable object
    df_table <- xtable(get(paste0(temp, "_tbl")), 
                       type = "latex", 
                       digits = 1)
    
    # Columns
    align(df_table) <- paste0("l", paste0(rep("c", length(loc_vec)), collapse = ""))
    
    # Header
    addtorow <- list()
    addtorow$pos <- list(0)
    addtorow$command <- paste0("Method & ", paste0(toupper(loc_vec), collapse = " & "), " \\\\\n")
    
    # Horizontal lines to include
    hline.after <- c(-1, 0, nrow(df_table))
    
    # Make Latex table
    print(file = paste0(fig_path, "table_qrd",  round(100*q_level, 2), "_", temp, ".txt"),
          x = df_table, 
          include.colnames = FALSE,
          add.to.row = addtorow,
          hline.after = hline.after,
          booktabs = TRUE,
          # only.contents = TRUE,
          format.args = list(big.mark = ",", decimal.mark = "."))
  }
}

#### PIT Histograms ####
# Number of bins in histogram
n_bins <- 12

# Limit for y-axis
y_lim <- 2

# List for matrices
df_plot <- data.frame(model = character(),
                      breaks = numeric(),
                      pit = numeric(),
                      stringsAsFactors = FALSE)

# For-Loop over methods
for(temp_meth in pp_vec){ 
  # Calculate histogram and read out values (see pit function)
  temp <- hist(pit_ls[[paste0(temp_meth, "_", loc_plot)]],
               breaks = (0:n_bins)/n_bins,
               plot = FALSE)
  
  # Get row of data frame
  i <- nrow(df_plot) + 1:(n_bins + 1)
  
  # Fill in data frame
  df_plot[i, "model"] <- temp_meth
  
  # Breaks
  df_plot[i, "breaks"] <- temp[["breaks"]]
  
  # Density
  df_plot[i, "pit"] <- c(0, temp[["density"]])
}

# Order locations
df_plot[["model"]] <- factor(df_plot[["model"]],
                             levels = pp_order,
                             labels = pp_labels)

# Start plot
pdf_plot <- 
  ggplot(df_plot, 
         aes(x = breaks, y = pit)) + 
  geom_rect(aes(xmin = breaks, xmax = lead(breaks),
                ymin = 0, ymax = lead(pit)),
            color = "black", fill = "black", alpha = 0.3) +
  facet_wrap("model",
             ncol = 4,
             labeller = label_parsed,
             scales = "fixed") +
  ylab("Density") +
  xlab("PIT") +
  theme_bw() +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("0", "0.25", "0.5", "0.75", "1")) +
  scale_y_continuous(breaks = c(0, 0.5, 1, 1.5, 2),
                     labels = c("0", "0.5", "1", "1.5", "2")) +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12)) + 
  geom_hline(aes(yintercept = 1), 
             linetype = "dashed") +
  geom_hline(aes(yintercept = y_lim), 
             linetype = "blank")

# Name of file
file_name <- paste0("pit_histograms_", loc_plot, ".pdf")

# Save as PDF
ggsave(filename = file_name,
       path = fig_path,
       plot = pdf_plot,
       width = 20,
       height = 5,
       scale = 0.7)

#### CORP reliability diagrams ####
# Quantile levels
q_levels <- c()

# Threshold levels
t_levels <- c(250)

# Query whether exceedance or non-exceedance is considered
log_exc <- FALSE

# Load data for calculation of quantile threshold
load(file = paste0(data_path, "data_", loc_plot, ".RData"))

# Thresholds to consider (based on quantiles and raw thresholds)
t_levels <- c(t_levels, quantile(x = data[["ghi"]],
                                 probs = q_levels,
                                 type = 8))

# Function that creates a panel of the reliability diagrams
fn_panel <- function(ls_rd, digits = 3){
  ### Input
  # ls_rd....List of reliability diagram objects (list of rd-objects)
  # digits...Digits for rounding of score decomposition (integer)
  ### Output
  # res...Panel plot (ggplot-object)
  ###
  
  #### Initiation ####
  # Read out models and quantiles of list
  model_vec <- pp_vec[which(sapply(pp_vec, function(x) any(grepl(x, names(ls_rd), fixed = TRUE))))]
  thr_vec <- t_levels[which(sapply(t_levels, function(x) 
    any(is.element(paste0(model_vec, round(x)), names(ls_rd)))))]
  
  # Columns of data frame
  col_vec <- c("model", "threshold" ,"x", "y", "x_rc" ,"lower", "upper", 
               "score", "mcb", "dsc", "unc")
  
  # Create data frame for panel
  df_rd <- as.data.frame(matrix(nrow = 0,
                                ncol = length(col_vec)))
  colnames(df_rd) <- col_vec
  
  #### Read out information ####
  # For-Loop over models and thresholds
  for(temp_model in model_vec){ for(thr in thr_vec){
    # Get element of list (if available)
    if(is.element(paste0(temp_model, round(thr)), names(ls_rd))){
      temp_rd <- ls_rd[[paste0(temp_model, round(thr))]][["p"]]
    }else{ next }
    
    # Rows of data frame to fill in
    i_df <- nrow(df_rd) + 1:nrow(temp_rd[["cases"]])
    
    # Read out information on plots
    df_rd[i_df, "model"] <- temp_model
    df_rd[i_df, "threshold"] <- round(thr)
    df_rd[i_df, "x"] <- temp_rd[["cases"]][["x"]]
    df_rd[i_df, "y"] <- temp_rd[["cases"]][["y"]]
    df_rd[i_df, "x_rc"] <- temp_rd[["cases"]][["CEP_pav"]]
    
    # Score decomposition based on rd-object
    temp_decomp <- summary(ls_rd[[paste0(temp_model, round(thr))]])
    
    # Read out scores
    df_rd[i_df, "score"] <- temp_decomp[["mean_score"]]
    df_rd[i_df, "mcb"] <- temp_decomp[["miscalibration"]]
    df_rd[i_df, "dsc"] <- temp_decomp[["discrimination"]]
    df_rd[i_df, "unc"] <- temp_decomp[["uncertainty"]]
    
    # Extra rows for lower and upper
    i_df <- nrow(df_rd) + 1:nrow(temp_rd[["regions"]])
    
    # Read out information of plots
    df_rd[i_df, "model"] <- temp_model
    df_rd[i_df, "threshold"] <- round(thr)
    df_rd[i_df, "x"] <- temp_rd[["regions"]][["x"]]
    df_rd[i_df, "lower"] <- temp_rd[["regions"]][["lower"]]
    df_rd[i_df, "upper"] <- temp_rd[["regions"]][["upper"]]
    
    # Read out scores
    df_rd[i_df, "score"] <- temp_decomp[["mean_score"]]
    df_rd[i_df, "mcb"] <- temp_decomp[["miscalibration"]]
    df_rd[i_df, "dsc"] <- temp_decomp[["discrimination"]]
    df_rd[i_df, "unc"] <- temp_decomp[["uncertainty"]]
  }}
  
  # Remove duplicates
  df_rd <- df_rd[!duplicated(df_rd),]
  
  # Group by (needed for score decomposition labeling)
  df_rd <- df_rd %>%
    group_by(model, threshold)
  
  # Order methods
  df_rd[["model"]] <- factor(df_rd[["model"]],
                             levels = pp_order,
                             labels = pp_labels)
  
  # Rename thresholds
  df_rd[["threshold"]] <- factor(df_rd[["threshold"]],
                                levels = t_levels,
                                labels = paste0(round(t_levels), " [W/m^2]"))
  
  #### Plot ####
  # Digits of score decomposition
  digits_f <- paste0("%.", digits, "f")
  
  # Scores of methods (Analogous to quantile reliability diagrams no S_bar)
  scores <- df_rd %>%
    distinct(across(score:unc)) %>%
    mutate(label = paste0(
      c("S ",  "MCB ", "DSC ", "UNC "),
      sprintf(digits_f, c(score, mcb, dsc, unc)),
      collapse = "\n"
    ))
  
  # Create layer with scores for plot
  score_layer <- list(
    geom_label(
      mapping = aes(x = -Inf, y = Inf, label = label),
      data = scores, size = 3.2, hjust = 0, vjust = 1, label.size = NA,
      alpha = 0, label.padding = unit(1, "lines"), parse = FALSE
    )
  )
  
  # Plotting together
  res <- ggplot(subset(df_rd, !is.na(x_rc)), aes(x, x_rc, group = model)) +
    geom_abline(intercept = 0, slope = 1, 
                colour = "grey70", size = 0.5) +
    geom_ribbon(data = subset(df_rd, !is.na(lower)), aes(ymin = lower, ymax = upper),
                color = "skyblue3", size = 0.6,
                fill = "skyblue3", alpha = 0.3) +
    geom_line(size = 1,
              color = "firebrick3") +
    xlab("Forecast Probability") +
    ylab("Conditional Event Probability") +
    scale_x_continuous(limits = c(0, 1),
                       breaks = c(0, 0.25, 0.5, 0.75, 1),
                       labels = c("0", "0.25", "0.5", "0.75", "1")) +
    scale_y_continuous(limits = c(0, 1),
                       breaks = c(0, 0.25, 0.5, 0.75, 1),
                       labels = c("0", "0.25", "0.5", "0.75", "1")) +
    theme_bw(base_size = 11) +
    coord_fixed() +
    theme(
      panel.grid.major = element_line(size = 0.05),
      panel.grid.minor = element_line(size = 0.05),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      strip.text.x = element_text(size = 12),
      strip.text.y = element_text(size = 12),
    ) +
    score_layer
  
  # Panel depending on number of models and quantiles
  if(length(model_vec) == 1){
    res <- res +
      facet_wrap("threshold",
                 ncol = 4,
                 scales = "fixed") +
      coord_fixed(ratio = 1)
  }else if(length(thr_vec) == 1){
    res <- res +
      facet_wrap("model",
                 ncol = 4,
                 scales = "fixed") +
      coord_fixed(ratio = 1)
  }else{
    res <- res +
      facet_grid(threshold ~ model) +
      coord_fixed(ratio = 1)
  }
  
  # Output
  return(res)
}

# Name of columns
col_vec <- c("date", "threshold", "value", "model", "truth")

# Make data frame
df <- as.data.frame(matrix(nrow = 0,
                           ncol = length(col_vec)))
colnames(df) <- col_vec

# For-Loop over methods
for(temp_meth in pp_vec){
  # Check whether local or global
  if(is.element(temp_meth, pp_local)){
    # Load data
    load(file = paste0(data_path, temp_meth, "_", loc_plot, ".RData"))
    
    # Create one data frame for each threshold
    fn_apply <- function(thr){
      # thr...Threshold
      
      # Make data frame
      res <- as.data.frame(matrix(nrow = nrow(data_te),
                                  ncol = length(col_vec)))
      colnames(res) <- col_vec
      
      # Read out date, location, model and truth (mod. for exc. later)
      res[,"date"] <- data_te[["Time"]]
      res[,"threshold"] <- thr
      res[,"model"] <- temp_meth
      res[,"truth"] <- as.numeric(data_te[["ghi"]] > thr)
      
      # Calculate exceedance probabilities (mod. for exc. later)
      if(temp_meth == "AnEn"){
        res[,"value"] <- rowMeans(pred[["f"]] > thr)
      }
      else if(grepl("idr", temp_meth, fixed = TRUE)){
        res[,"value"] <- 1 - as.vector(isodistrreg::cdf(predictions = pred[["pred_idr"]],
                                                        thresholds = thr))
      }
      
      # Output
      return(res)
    }
    
    # Apply on quantile levels
    df <- rbind(df, bind_rows(lapply(t_levels, fn_apply)))
  }
  else if(is.element(temp_meth, pp_global)){
    # Get name of file
    if(((temp_meth == "drn") & drn_bias) |
       ((temp_meth == "bqn") & bqn_bias)){ temp_name <- paste0(temp_meth, "_bias") 
    }else{ temp_name <- temp_meth }
    
    # Load data
    load(file = paste0(data_path, temp_name, ".RData"))
    
    # Subset of data
    i_loc <- which(data_te[["location"]] == loc_plot)
    
    # Subset of data
    data_sub <- subset(data_te, location == loc_plot)
    
    # Create one data frame for each threshold
    fn_apply <- function(thr){
      # thr...Threshold
      
      # Make data frame
      res <- as.data.frame(matrix(nrow = nrow(data_sub),
                                  ncol = length(col_vec)))
      colnames(res) <- col_vec
      
      # Read out date, location, model and truth (mod. for exc. later)
      res[,"date"] <- data_sub[["Time"]]
      res[,"threshold"] <- thr
      res[,"model"] <- temp_meth
      res[,"truth"] <- as.numeric(data_sub[["ghi"]] > thr)
      
      # Calculate exceedance probabilities (mod. for exc. later)
      if(temp_meth == "drn"){
        res[,"value"] <- 1 - crch::ptnorm(q = thr, 
                                          mean = pred[["f"]][i_loc, 1], 
                                          sd = pred[["f"]][i_loc, 2],
                                          left = 0)
      }
      else if(temp_meth == "bqn"){
        res[,"value"] <- sapply(i_loc, function(i){
          mean(pred[["f"]][i,] > thr) })
      }
      
      # Output
      return(res)
    }
    
    # Apply on quantile levels
    df <- rbind(df, bind_rows(lapply(t_levels, fn_apply)))
    rm(data_sub)
  }
  rm(pred, data_te)
}

# Modify for non-exceedance if needed
if(!log_exc){
  df[["value"]] <- 1 - df[["value"]]
  df[["truth"]] <- 1 - df[["truth"]]
}

# For-Loop over thresholds
for(thr in t_levels){
  # List for reliability diagrams
  ls_rd <- list()
  
  # For-Loop over methods
  for(temp_meth in pp_order){
    # Subset
    df_sub <- subset(df, (model == temp_meth) & 
                       (threshold == thr))
    
    # Make reliability diagram of predictions
    ls_rd[[paste0(temp_meth, round(thr))]] <- 
      reliabilitydiag(p = df_sub[["value"]],
                      y = df_sub[["truth"]],
                      region.level = 0.9,
                      region.position = "diagonal",
                      n.boot = 1000)
  }
  
  # Plot together
  fig_plot <- fn_panel(ls_rd)
  
  # Name dependent on location
  if(log_exc){ file_name <- paste0(fig_path, "threshold_exceedance_corp_t", 
                                   round(thr), "_", loc_plot, ".pdf") }
  else{ file_name <- paste0(fig_path, "threshold_non-exceedance_corp_t", 
                                   round(thr), "_", loc_plot, ".pdf") }
  
  # Save
  ggplot2::ggsave(filename = file_name,
                  plot = fig_plot,
                  width = 14,
                  height = 5*ceiling(length(pp_vec)/4),
                  scale = 0.7
  )
}

#### Quantile reliability diagram ####
# Quantile level
q_levels <- c(1, 11)/12

# Number of resamples
n_resamples <- 999

# Include points
log_points <- TRUE

# Load data frame
load(file = paste0(data_path, "data_quantile_rd.RData"))

# Subset of data
df_sub <- subset(df, is.element(quantile, q_levels) & 
                   is.element(model, pp_vec) & 
                   (location == toupper(loc_plot)))

# Name of file for quantile reliability diagram data
file_name <- paste0(data_path, "qreldiag_", n_resamples, "_", paste0(pp_vec, collapse = "_"), "_q", 
                    paste0(round(100*q_levels, 2), collapse = "_"), "_", loc_plot, ".RData")

# Check if file exists
if(file.exists(file_name)){ load(file = file_name)
}else{
  # Data frame for plotting
  df_reldiag <- df_sub %>%
    group_by(model, quantile) %>%
    summarize(reldiag(value, truth, alpha = unique(quantile),
                      n_resamples = n_resamples, digits = 2, ties = "primary"),
              .groups = "keep") %>%
    mutate(across(c(x_rc, lower, upper), ~ pmax(., 0))) # set negative values to zero
  
  # Save data
  save(file = file_name,
       list = c("df_reldiag"))
}

# Get maximum of REST2 and GHI?
load(file = paste0(data_path, "data_bon.RData"))
temp_max <- max(subset(data, year(Time) == 2020)[, c("ghi", "rest2")])
rm(data)

# Cut consistency bands
for(temp_bd in c("lower", "upper")){
  df_reldiag[[temp_bd]] <- pmax(0, df_reldiag[[temp_bd]])
  df_reldiag[[temp_bd]] <- pmin(df_reldiag[[temp_bd]], temp_max)
}

# Order methods
df_reldiag[["model"]] <- factor(df_reldiag[["model"]],
                                levels = pp_order,
                                labels = pp_labels)

# Plot
pdf_plot <-
  plot_reldiag(df_reldiag, mcb_decomp = FALSE, pval = FALSE, plot_pts = log_points,
               square = TRUE, my_alpha = 0.1, score_digits = 1) +
  facet_grid(quantile ~ model,
             scales = "fixed") +
  ylab(NULL) +
  coord_fixed(ratio = 1)

# Name of file
temp_name <- paste0(paste0("quantile_reldiag", n_resamples, "_q", 
                           paste0(round(100*q_levels, 2), collapse = "_"), "_", loc_plot ,".pdf"))

# Save as PDF
ggsave(filename = temp_name,
       path = fig_path,
       plot = pdf_plot,
       width = 14,
       height = 8,
       scale = 0.7)
