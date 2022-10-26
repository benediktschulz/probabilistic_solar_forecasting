## Functions for evaluation of postprocessed forecasts ##

#### Import ####
# Import basic functions
source(paste0(getwd(), "/functions_basic.R"))

#### Coverage ####
# Calculate coverage of a central prediction interval
fn_cover <- function(x, alpha = 0.1){
  ###-----------------------------------------------------------------------------
  ###Input
  #x.......PIT values (n vector)
  #alpha...Significance level (probability)
  #........Default: 0.1 -> 10% -> 90% prediction interval
  ###-----------------------------------------------------------------------------
  ###Output
  #res...Coverage in percentage
  ###-----------------------------------------------------------------------------
  
  #### Coverage calculation ####
  res <- mean((alpha/2 <= x) & (x <= (1 - alpha/2)))

  # Output as percentage
  return(100*res)
}

#### Brier score ####
# Brier score for given distribution or ensemble
brier_score <- function(f, y, t = 0, distr = "ens", t_distr = 0){
  ###-----------------------------------------------------------------------------
  ###Input
  #f.........distr == par. distr.: Matrix with location and scale of forecast distribution (n x n_par matrix)
  #..........distr == "ens": Ensemble forecasts (n x n_ens matrix)
  #..........distr == "p" (or elsewise): Probability forecasts of exceeding t (n vector)
  #y.........Observations (n vector)
  #t.........Brier Score Threshold (non-negative scalar)
  #..........Default: 0
  #distr.....Forecast type (specific distribution or ensemble) (string)
  #..........Default: "ens" -> Ensemble
  #t_distr...Threshold for censored or truncated distribution (Scalar)
  #..........Default: 0
  ###-----------------------------------------------------------------------------
  ###Output
  #res...Brier Scores of n forecasts (n vector)
  ###-----------------------------------------------------------------------------
  
  #### Calculation ####
  ## Brier Score w.r.t. event of exceeding threshold t
  ## p_t = 1 - F(t). BS_t (F, y) = (p_t - 1(y > t))^2 = (F(t) - 1(y <= t))^2
  
  # Calculate F(t) depending on distribution
  if(distr == "ens"){ 
    if(is.vector(f)){ f <- mean(f <= t) }
    else{ f <- rowMeans(f <= t) }
  }
  # Truncated logistic
  else if(distr == "tlogis"){ f <- (t > t_distr)*crch::ptlogis(q = t, 
                                                               location = f[,1], 
                                                               scale = f[,2],
                                                               left = t_distr) }
  # Truncated normal
  else if(distr == "tnorm"){ f <- (t > t_distr)*crch::ptnorm(q = t, 
                                                             mean = f[,1], 
                                                             sd = f[,2],
                                                             left = t_distr) }
  
  # Calculate Brier Score
  res <- (f - (y <= t))^2
  
  # Return score
  return(res)
}

#### Pinball loss ####
# Quantile score / pinball loss for given distribution or ensemble
pinball_loss <- function(f, y, alpha = 0.95, distr = "ens", t_distr = 0){
  ###-----------------------------------------------------------------------------
  ###Input
  #f.........distr == par. distr.: Matrix with location and scale of forecast distribution (n x n_par matrix)
  #..........distr == "ens": Ensemble forecasts (n x n_ens matrix)
  #..........distr == "q" (or elsewise): Quantile forecasts at level alpha (n vector)
  #y.........Observations (n vector)
  #alpha.....Quantile level (probability)
  #..........Default: 0.95 -> 95%
  #distr.....Forecast type (specific distribution or ensemble) (string)
  #..........Default: "ens" -> Ensemble
  #t_distr...Threshold for censored or truncated distribution (Scalar)
  #..........Default: 0
  ###-----------------------------------------------------------------------------
  ###Output
  #res...Quantile scores / pinball losses of n forecasts (n vector)
  ###-----------------------------------------------------------------------------
  
  #### Calculation ####
  ## Pinball loss of quantile forecast q_f at level alpha
  ## q_f = Q(alpha). PL (q_f, y) = (q_f - y)*(1(y <= q_f) - alpha)
  
  # Calculate F(t) depending on distribution
  if(distr == "ens"){ 
    if(is.vector(f)){ f <- quantile(x = f,
                                    probs = alpha,
                                    type = 8) }
    else{ f <- sapply(1:nrow(f), function(x) quantile(x = f[x,],
                                                      probs = alpha,
                                                      type = 8)) }
  }
  # Truncated logistic
  else if(distr == "tlogis"){ f <- crch::qtlogis(p = alpha, 
                                                 location = f[,1], 
                                                 scale = f[,2],
                                                 left = t_distr) }
  # Truncated normal
  else if(distr == "tnorm"){ f <- crch::qtnorm(p = alpha, 
                                               mean = f[,1], 
                                               sd = f[,2],
                                               left = t_distr) }
  
  # Calculate quantile score / pinball loss
  res <- (f - y)*((y <= f) - alpha)
  
  # Return score
  return(res)
}

#### Interval score ####
# Interval score of central prediction interval based on upper and lower boundaries
interval_score <- function(l, u, y, alpha = 0.1){
  ###-----------------------------------------------------------------------------
  ###Input
  #l.......Lower boundaries of prediction interval (n vector)
  #u.......Upper boundaries of prediction interval (n vector)
  #y.......Observations (n vector)
  #alpha...1 - Level of prediction interval (probability)
  #........Default: 0.1 -> 90% prediction interval
  ###-----------------------------------------------------------------------------
  ###Output
  #res...Interval scores of n forecasts (n vector)
  ###-----------------------------------------------------------------------------
  
  #### Calculation ####
  # Calculate interval score
  res <- (u - l) + 2/alpha*(l - y)*(y < l) + 2/alpha*(y - u)*(y > u)
  
  # Return score
  return(res)
}

#### BQN: Bernstein Quantile function ####
# Function that calculates quantiles for given coefficients
bern_quants <- function(alpha, q_levels){
  ###-----------------------------------------------------------------------------
  ###Input
  #alpha......Coefficients of Bernstein Basis (n x (p_degree + 1) matrix)
  #q_levels...Quantile levels (n_q vector)
  ###-----------------------------------------------------------------------------
  ###Output
  #res...Quantile forecasts for given coefficients (n x n_q matrix)
  ###-----------------------------------------------------------------------------
  
  #### Initiation ####
  # Get degree of polynomials from coefficients
  if(is.vector(alpha)){ p_degree <- length(alpha) - 1 }
  else{ p_degree <- ncol(alpha) - 1 }
  
  #### Calculation ####
  # Calculate quantiles (sum of coefficients times basis polynomials)
  if(length(q_levels) == 1){ res <- alpha %*% sapply(0:p_degree, dbinom, size = p_degree, prob = q_levels) }
  else{ res <- alpha %*% t(sapply(0:p_degree, dbinom, size = p_degree, prob = q_levels)) }

  # Return quantiles
  return(res)
}

#### Evaluation of ensemble forecasts ####
# Function to calculate evaluation measures of scores
fn_scores_ens <- function(ens, y, alpha = 0.1, skip_evals = NULL, scores_ens = TRUE){
  ###-----------------------------------------------------------------------------
  ###Input
  #ens..........Ensemble data for prediction (n x n_ens matrix)
  #y............Observations for prediction (n vector)
  #alpha........1 - Level of prediction interval (probability)
  #.............Default: 0.1 -> 90% prediction interval
  #skip_evals...Skip the given evaluation measures (string vector)
  #.............Default: NULL -> Calculate all
  #scores_ens...Should scores of ensemble forecasts be calculated? (logical)
  #.............Default: TRUE
  ###-----------------------------------------------------------------------------
  ###Output
  #...scores_ens...Data frames containing (n x 6 data frame):
  #......rank......Ranks of observations in ensemble forecasts (n vector)
  #......crps......CRPS of ensemble forecasts (n vector)
  #......logs......Log-Score of ensemble forecasts (n vector)
  #......lgt.......Ensemble range (n vector)
  #......e_md......Bias of median forecast (n vector)
  #......e_me......Bias of mean forecast (n vector)
  ###-----------------------------------------------------------------------------
  
  #### Initiation ####
  # Load packages
  library(scoringRules)
  
  # Calculate only if scores_ens is TRUE
  if(!scores_ens){ return(FALSE) }
  
  # Check if vector is given
  if(is.vector(ens)){ ens <- matrix(data = ens,
                                    nrow = 1) }
  
  # Get number of ensembles
  n <- nrow(ens)
  
  # Get ensemble size
  n_ens <- ncol(ens)
  
  # Make data frame
  scores_ens <- data.frame(rank = numeric(length = n),
                           crps = numeric(length = n),
                           logs = numeric(length = n),
                           lgt = numeric(length = n),
                           e_me = numeric(length = n),
                           e_md = numeric(length = n))
  
  #### Calculation ####
  # Calculate observation ranks
  if(is.element("rank", colnames(scores_ens))){
    scores_ens[["rank"]] <- apply(cbind(y, ens), 1, function(x){ rank(x, ties = "random")[1] }) }
  
  # Calculate CRPS of raw ensemble
  if(is.element("crps", colnames(scores_ens))){
    scores_ens[["crps"]] <- crps_sample(y = y, 
                                        dat = ens) }
  
  # Calculate Log-Score of raw ensemble
  if(is.element("logs", colnames(scores_ens))){
    scores_ens[["logs"]] <- logs_sample(y = y, 
                                        dat = ens) }
  
  # Calculate (1-alpha) prediction interval
  if(is.element("lgt", colnames(scores_ens))){
    # Corresponding quantiles (alpha/2 and 1-alpha/2) are included
    if( ((((n_ens + 1)*alpha/2) %% 1) == 0) & 
       ((((n_ens + 1)*(1-alpha/2)) %% 1) == 0) ){ 
      # Indices of corresponding quantiles
      i_lgt <- (n_ens + 1)*c(alpha/2, 1-alpha/2)
      
      # Get quantiles
      q_lgt <- t(apply(ens, 1, sort))[,i_lgt]
      
      # Transform if vector
      if(n == 1){ q_lgt <- matrix(data =  q_lgt,
                                  nrow = 1) }
      
      # Calculate corresponding range
      scores_ens[["lgt"]] <- apply(t(apply(q_lgt, 1, range)), 1, diff) 
    }
    # Quantiles are not included: Calculate corresponding via quantile function
    else{ 
      # Choose type 8, as suggested in ?quantile (and bias observed for linearly pooled aggregation)
      scores_ens[["lgt"]] <- apply(ens, 1, function(x) 
        diff(quantile(x = x,
                      probs = c(alpha/2, 1-alpha/2),
                      type = 8)) ) 
    }
  }
  
  # Calculate bias of median forecast
  if(is.element("e_md", colnames(scores_ens))){
    scores_ens[["e_md"]] <- apply(ens, 1, median) - y }
  
  # Calculate bias of mean forecast
  if(is.element("e_me", colnames(scores_ens))){
    scores_ens[["e_me"]] <- rowMeans(ens) - y }
  
  #### Output ####
  # Skip evaluation measures
  scores_ens <- as.data.frame(scores_ens[,!is.element(colnames(scores_ens), skip_evals), drop = FALSE])
  
  # Return output
  return(scores_ens)
}

#### Evaluation of parametric distributional forecasts ####
# Function for prediction based on the distributional parameters #
fn_scores_distr <- function(f, y, distr = "tlogis", alpha = 0.1,
                            skip_evals = NULL){
  ###-----------------------------------------------------------------------------
  ###Input
  #f............Parameters of forecast distribution (n x n_par matrix)
  #y............Observations (n vector)
  #distr........Parametric distribution ("tlogis", "tnorm", "norm")
  #.............Default: (zero-)truncated logistic
  #alpha........1 - Level of prediction interval (probability)
  #.............Default: 0.1 -> 90% prediction interval
  #skip_evals...Skip the following evaluation measures (string vector)
  #.............Default: NULL -> Calculate all
  ###-----------------------------------------------------------------------------
  ###Output
  #res...List containing:
  #......scores_pp...Data frames containing (n x 6 data frame):
  #.........pit.........PIT values of distributional forecasts (n vector)
  #.........crps........CRPS of forecasts (n vector)
  #.........logs........Log-Score of forecasts (n vector)
  #.........lgt.........Length of prediction interval (n vector)
  #.........e_md........Bias of median forecast (n vector)
  #.........e_me........Bias of mean forecast (n vector)
  ###-----------------------------------------------------------------------------
  
  #### Initiation ####
  # Load packages
  library(scoringRules)
  
  # Input check
  if(!is.element(distr, c("tlogis", "tnorm", "norm"))){ 
    print("Chosen 'distr' not available. Choose 'tlogis', 'tnorm' or 'norm.") }
  if(is.element(distr, c("tlogis", "tnorm", "norm")) & any(f[,2] < 0)){ print("Non-positive scale forecast!") }
  
  #### Data preparation ####
  # Number of predictions
  n <- nrow(f)
  
  # Make data frame
  scores_pp <- data.frame(pit = numeric(length = n),
                          crps = numeric(length = n),
                          logs = numeric(length = n),
                          lgt = numeric(length = n),
                          e_me = numeric(length = n),
                          e_md = numeric(length = n))
  
  #### Prediction and score calculation ####
  # Forecasts depending on distribution
  if(distr == "tlogis"){ # truncated logistic
    # Calculate PIT values
    if(is.element("pit", colnames(scores_pp))){
      scores_pp[["pit"]] <- crch::ptlogis(q = y, 
                                          location = f[,1], 
                                          scale = f[,2],
                                          left = 0) }
    
    # Calculate CRPS of forecasts
    if(is.element("crps", colnames(scores_pp))){
      scores_pp[["crps"]] <- crps_tlogis(y = y, 
                                         location = f[,1], 
                                         scale = f[,2],
                                         lower = 0) }
    
    # Calculate Log-Score of forecasts
    if(is.element("logs", colnames(scores_pp))){
      scores_pp[["logs"]] <- logs_tlogis(y = y, 
                                         location = f[,1], 
                                         scale = f[,2],
                                         lower = 0) }
    
    # Calculate length of (1 - alpha) % prediction interval
    if(is.element("lgt", colnames(scores_pp))){
      scores_pp[["lgt"]] <- crch::qtlogis(p = (1 - alpha/2), 
                                          location = f[,1], 
                                          scale = f[,2],
                                          left = 0) - crch::qtlogis(p = alpha/2, 
                                                                    location = f[,1], 
                                                                    scale = f[,2],
                                                                    left = 0) }
    
    # Calculate bias of median forecast
    if(is.element("e_md", colnames(scores_pp))){
      scores_pp[["e_md"]] <- crch::qtlogis(p = 0.5, 
                                           location = f[,1], 
                                           scale = f[,2],
                                           left = 0) - y }
    # scores_pp[["e_md"]] <- (f[,1] + f[,2]*log(1 + 2*exp(- f[,1]/f[,2]))) - y
    
    # Calculate bias of mean forecast
    if(is.element("e_me", colnames(scores_pp))){
      scores_pp[["e_me"]] <- (f[,1] - f[,2]*log(1 - plogis(- f[,1]/f[,2])))/(1 - plogis(- f[,1]/f[,2])) - y }
  }
  else if(distr == "tnorm"){ # truncated normal
    # Calculate PIT values
    if(is.element("pit", colnames(scores_pp))){
      scores_pp[["pit"]] <- crch::ptnorm(q = y, 
                                         mean = f[,1], 
                                         sd = f[,2],
                                         left = 0) }
    
    # Calculate CRPS of forecasts
    if(is.element("crps", colnames(scores_pp))){
      scores_pp[["crps"]] <- crps_tnorm(y = y, 
                                        location = f[,1], 
                                        scale = f[,2],
                                        lower = 0) }
    
    # Calculate Log-Score of forecasts
    if(is.element("logs", colnames(scores_pp))){
      scores_pp[["logs"]] <- logs_tnorm(y = y, 
                                        location = f[,1], 
                                        scale = f[,2],
                                        lower = 0) }
    
    # Calculate length of (1 - alpha) % prediction interval
    if(is.element("lgt", colnames(scores_pp))){
      scores_pp[["lgt"]] <- crch::qtnorm(p = (1 - alpha/2), 
                                         mean = f[,1], 
                                         sd = f[,2],
                                         left = 0) - crch::qtnorm(p = alpha/2, 
                                                                  mean = f[,1], 
                                                                  sd = f[,2],
                                                                  left = 0) }
    
    # Calculate bias of median forecast
    if(is.element("e_md", colnames(scores_pp))){
      scores_pp[["e_md"]] <- crch::qtnorm(p = 0.5, 
                                          mean = f[,1], 
                                          sd = f[,2],
                                          left = 0) - y }
    
    # # Calculate bias of mean forecast
    # scores_pp[["e_me"]] <- #TODO
  }
  else if(distr == "norm"){ # normal
    # Calculate PIT values
    if(is.element("pit", colnames(scores_pp))){
      scores_pp[["pit"]] <- pnorm(q = y, 
                                  mean = f[,1], 
                                  sd = f[,2]) }
    
    # Calculate CRPS of forecasts
    if(is.element("crps", colnames(scores_pp))){
      scores_pp[["crps"]] <- crps_norm(y = y, 
                                       location = f[,1], 
                                       scale = f[,2]) }
    
    # Calculate Log-Score of forecasts
    if(is.element("logs", colnames(scores_pp))){
      scores_pp[["logs"]] <- logs_norm(y = y, 
                                       location = f[,1], 
                                       scale = f[,2]) }
    
    # Calculate length of (1 - alpha) % prediction interval
    if(is.element("lgt", colnames(scores_pp))){
      scores_pp[["lgt"]] <- qnorm(p = (1 - alpha/2), 
                                  mean = f[,1], 
                                  sd = f[,2]) - qnorm(p = alpha/2, 
                                                      mean = f[,1], 
                                                      sd = f[,2]) }
    
    # Calculate bias of median forecast
    if(is.element("e_md", colnames(scores_pp))){
      scores_pp[["e_md"]] <- qnorm(p = 0.5, 
                                   mean = f[,1], 
                                   sd = f[,2]) - y }
    
    # Calculate bias of mean forecast
    scores_pp[["e_me"]] <- f[,1] - y
  }
  
  #### Output ####
  # Skip evaluation measures
  scores_pp <- as.data.frame(scores_pp[,!is.element(colnames(scores_pp), skip_evals), drop = FALSE])
  
  # Return
  return(scores_pp)
}