## Functions for Postprocessing ##

#### Import ####
# Import basic functions
source(paste0(getwd(), "/functions_basic.R"))
source(paste0(getwd(), "/functions_data.R"))
source(paste0(getwd(), "/functions_eval.R"))

#### IDR ####
# Function for estimating IDR (with subagging) #
idr_pp <- function(train, X, obs_var = "ghi", pred_vars = c("ssrd"), 
                   q_levels = NULL, idr_ls = list(), scores_pp = TRUE){
  ###-----------------------------------------------------------------------------
  ###Input
  #train.........Training data including predictors and obs. (n_train x (n_preds + 1) data.frame)
  #X.............Test data (-"-)
  #pred_vars.....Predictors used for IDR (n_preds vector of strings)
  #..............Default: c("ssrd") -> Only ssrd
  #obs_var.......Name of observed variable (string)
  #..............Default: "ghi"
  #q_levels......Quantile levels used for output and evaluation (probability vector)
  #..............Default: NULL -> Steps of 0.01
  #idr_ls........List that may contain the following variables:
  #...groups.....Groups of predictors corresponding to a partial order (named vector)
  #..............Default: NULL -> All predictors in one group
  #...orders.....Define partial orders of given groups (named vector)
  #..............Default: NULL -> All increasing convex order (exchangeable ensemble)
  #...n_sbg......Number of subsamples for subagging (integer)
  #..............Default: 0 -> No subagging
  #...n_sample...Size of subsamples in subagging (integer)
  #..............Default: Half of training set, at most 1,000
  #...max_iter...OSQP parameter: Number of maximum iterations (integer)
  #..............Default: 1000L
  #...eps_abs....OSQP parameter: Absolute convergence tolerance
  #..............Default: 1e-3
  #...eps_rel....OSQP parameter: Relative convergence tolerance
  #..............Default: 1e-3
  #...console....Show output of IDR function (logical)
  #..............Default: FALSE (Use FALSE and not 0 !!)
  #scores_pp.....Should evaluation measures be calculated (logical)
  #..............Default: TRUE
  ###-----------------------------------------------------------------------------
  ###Output
  #res...List containing:
  #......f..............IDR quantile forecasts (n x n_q matrix)
  #......q_levels.......Quantile levels of forecasts in f (n_q probability vector)
  #......idr_ls.........Hyperparameters (list)
  #......runtime_est....Estimation time (no sbg), incl. prediction (sbg) (numeric)
  #......runtime_pred...Prediction time (no sbg) (numeric)
  #......n_train........Number of training samples (integer)
  #......n_test.........Number of test samples (integer)
  #......scores_pp......Data frame containing (n x 6 data frame):
  #.........pit.........PIT of IDR forecasts (n vector)
  #.........crps........CRPS of IDR forecasts (n vector)
  #.........logs........Log-Score of IDR forecasts (n vector)
  #.........lgt.........Length of IDR prediction interval (n vector)
  #.........e_md........Bias of median forecast (n vector)
  #.........e_me........Bias of mean forecast (n vector)
  ###-----------------------------------------------------------------------------
  
  #### Initiation ####
  # Load packages
  library(isodistrreg)
  library(scoringRules)
  
  # Relevant variables for training and testing
  train_vars <- c(obs_var, pred_vars)
  test_vars <- pred_vars
  
  # Observations for scores
  if(scores_pp){ test_vars <- c(test_vars, obs_var) }
  
  # Input check
  if(any(!is.element(train_vars, names(train)))){ 
    print("IDR: Training data does not include relevant variables.") }
  if(any(!is.element(test_vars, names(X)))){ 
    print("IDR: Test data does not include all of the relevant variables.") }
  
  # Cut test/training data to relevant variables
  train <- train[,train_vars]
  X <- X[,test_vars]
  
  # Input check
  if(any(is.na(train))){
    print("IDR: Training data includes missing values! Missing values are left out!")
    train <- na.omit(train)
  }
  if(any(is.na(X))){
    print("IDR: Test data includes missing values! Missing values are left out!")
    X <- na.omit(X)
  }
  
  # If not given use equidistant quantiles (steps of 0.01, incl. median)
  if(is.null(q_levels)){ q_levels <- seq(from = 1/100,
                                         to = 99/100,
                                         by = 1/100) }
  
  #### Hyperparameter ####
  # Hyperparameters and their default values
  hpar_ls <- list(console = FALSE,
                  groups = rep(1, length(pred_vars)),
                  orders = c("sd" = 1),
                  n_sbg = 0,
                  n_sample = min(1000, floor(nrow(train)/2)),
                  max_iter = 1000L,
                  eps_abs = 1e-3,
                  eps_rel = 1e-3)
  
  # Update hyperparameters
  idr_ls <- update_hpar(hpar_ls = hpar_ls,
                        in_ls = idr_ls)
  
  # Set names for groups
  idr_ls$groups <- setNames(idr_ls$groups, pred_vars)

  # Define optimization settings
  pars_osqp <- list(verbose = FALSE,
                    eps_abs = idr_ls$eps_abs,
                    eps_rel = idr_ls$eps_rel,
                    max_iter = idr_ls$max_iter)
  
  #### Data preparation ####
  # Number of predictions
  n <- nrow(X)
  
  #### Estimation/Prediction ####
  # Differentiate subagging
  if(idr_ls$n_sbg == 0){ 
    # Take time
    start_tm <- Sys.time()
    
    # Estimate IDR fit
    est <- idr(y = train[[obs_var]],
               X = train[,pred_vars],
               groups = idr_ls$groups,
               orders = idr_ls$orders,
               progress = idr_ls$console,
               pars = pars_osqp)
    
    # Take time
    end_tm <- Sys.time()
    
    # Time needed
    runtime_est <- as.numeric(difftime(end_tm, start_tm, units = "mins"))
    
    # Take time
    start_tm <- Sys.time()
    
    # Calculate prediction object
    pred_idr <- predict(object = est,
                        data = X)
    
    # Take time
    end_tm <- Sys.time()
    
    # Time needed
    runtime_pred <- as.numeric(difftime(end_tm, start_tm, units = "mins"))
  }
  else{ 
    # Take time
    start_tm <- Sys.time()
    
    # Estimate and predict via subagging
    pred_idr <- idrbag(y = train[[obs_var]],
                       X = train[,pred_vars],
                       groups = idr_ls$groups, 
                       orders = idr_ls$orders,
                       pars = pars_osqp, 
                       progress = idr_ls$console, 
                       newdata = X[,pred_vars], 
                       b = idr_ls$n_sbg, 
                       p = idr_ls$n_sample/nrow(train),
                       replace = (nrow(train) < idr_ls$n_sample*idr_ls$n_sbg)) 
    
    # Take time
    end_tm <- Sys.time()
    
    # Time needed
    runtime_est <- as.numeric(difftime(end_tm, start_tm, units = "mins"))
    
    # No specific prediction runtime
    runtime_pred <- NA
  }
  
  #### Quantiles and Evaluation ####
  # Predict quantiles
  q <- qpred(predictions = pred_idr,
             quantiles = q_levels)
  
  # Level of prediction interval
  alpha <- 0.1
  
  # Initialize evaluation measure of IDR forecasts
  if(scores_pp){ 
    scores_pp <- data.frame(pit = numeric(length = n),
                            crps = numeric(length = n),
                            logs = numeric(length = n),
                            lgt = numeric(length = n),
                            e_me = numeric(length = n),
                            e_md = numeric(length = n))
    
    
    # Calculate PIT values
    scores_pp[["pit"]] <- pit(predictions = pred_idr,
                              y = X[[obs_var]])
    
    
    # Calculate CRPS of IDR forecasts
    scores_pp[["crps"]] <- isodistrreg::crps(predictions = pred_idr,
                                             y = X[[obs_var]])
    
    # Calculate bias of median forecast
    if(any(q_levels == 0.5)){ scores_pp[["e_md"]] <- q[,which(q_levels == 0.5)] - X[[obs_var]] }
    else{ scores_pp[["e_md"]] <- qpred(predictions = pred_idr,
                                       quantiles = 0.5) - X[[obs_var]] }
    
    # Calculate length of 90% IDR prediction interval (based on given quantiles)
    if(any(q_levels == alpha/2) & any(q_levels == (1 - alpha/2))){ 
      scores_pp[["lgt"]] <- q[,which(q_levels == (1 - alpha/2))] - q[,which(q_levels == alpha/2)] }
    # Calculate quantiles if not given
    else{ scores_pp[["lgt"]] <- apply(qpred(predictions = pred_idr,
                                            quantiles = c(alpha/2, (1 - alpha/2))), 1, diff) }
    
    # Calculate Log-Score and bias of mean forecast as for ensemble
    temp_scores <- fn_scores_ens(ens = q,
                                 y = X[[obs_var]],
                                 skip_evals = c("rank", "crps", "e_md", "lgt"),
                                 scores_ens = TRUE)
    
    # Read out Log-Score and bias of mean forecast
    scores_pp[["logs"]] <- temp_scores[["logs"]]
    scores_pp[["e_me"]] <- temp_scores[["e_me"]]
  }
  
  #### Output ####
  return(list(f = q, 
              q_levels = q_levels,
              pred_idr = pred_idr,
              idr_ls = idr_ls,
              scores_pp = scores_pp, 
              pred_vars = pred_vars,
              obs_var = obs_var,
              n_preds = length(pred_vars),
              n_train = nrow(train),
              n_test = nrow(X),
              runtime_est = runtime_est,
              runtime_pred = runtime_pred))
}

#### DRN ####
## Distribution: (zero-)truncated logistic/normal distribution
## Estimation: CRPS
## Network ensemble aggregation: Parameter averaging (Quantile averaging, Vincentization)
## Comment: Station embedding

# Function for pp including estimation and prediction
drn_pp <- function(train, X, i_valid = NULL, loc_id_vec = NULL,
                   pred_vars = c("ssrd", "location"), obs_var = "ghi",
                   nn_ls = list(), n_cores = NULL, model_out = FALSE,
                   scores_pp = TRUE){
  ###-----------------------------------------------------------------------------
  ###Input
  #train...........Training data including predictors and obs. (n_train x (n_preds + 1) data.frame)
  #X...............Test data including predictors (and obs.) (-"-)
  #i_valid.........Indices of validation data (n_valid vector)
  #................Default: NULL -> Use fraction of training data
  #loc_id_vec......IDs of locations (string vector)
  #................Default: NULL -> Only needed if locations are included
  #pred_vars.......Predictors used for NN (vector of strings)
  #................Default: c("ssrd", "location")
  #obs_var.........Name of observed variable (string)
  #................Default: "ghi"
  #nn_ls...........List that may contain the following variables:
  #...n_sim........Number of NN to estimate and average (positive integer)
  #................Default: 10
  #...distr........Forecast distribution ("tlogis", "tnorm")
  #................Default: "tnorm"
  #...t_scale......Threshold for scale parameter truncation (positive scalar)
  #................Default: 1e-4
  #...lr_adam......Learning rate in Adam-Optimizer (scalar)
  #................Default: 5e-4
  #...n_epochs.....Number of epochs (integer)
  #................Default: 150
  #...n_patience...Patience in callbacks (integer)
  #................Default: 10
  #...n_batch......Size of batches is equal to n_train/n_batch (input parameter batch_size)
  #................Default: 64
  #...emb_dim......Dimension of station embedding (scalar)
  #................Default: 10
  #...layers.......Vector that defines the nodes in the hidden layers (n_layers vector)
  #................Default: c(64, 32) -> 2 layers: 64 in first, 32 in second
  #...actv.........Activation function of non-output layers (string)
  #................Default: "softplus"
  #...nn_verbose...Query, if network output should be shown (logical)
  #................Default: 0
  #n_cores.........Number of cores used in keras (integer)
  #................Default: NULL -> Use one less than available
  #model_out.......Query, if models should be in output (logical)
  #................Default: FALSE -> No model output
  #scores_pp.......Should evaluation measures be calculated (logical)
  #................Default: TRUE
  ###-----------------------------------------------------------------------------
  ###Output
  #res...List containing:
  #......f..............Distributional forecasts (n x 2 matrix)
  #......nn_ls..........Hyperparameters (list)
  #......model_ls.......Network models of the ensemble (list)
  #......pred_vars......Predictors (string vector)
  #......n_preds........Number of predictors (integer)
  #......n_train........Number of training samples (integer)
  #......n_valid........Number of validation samples (integer)
  #......n_test.........Number of test samples (integer)
  #......runtime_est....Estimation time (numeric)
  #......runtime_pred...Prediction time (numeric)
  #......center.........Centering parameter of training data (numeric)
  #......scale..........Scaling parameter of training data (numeric)
  #......scores_pp......Data frames containing (n x 6 data frame):
  #.........pit.........PIT of DRN forecasts (n vector)
  #.........crps........CRPS of DRN forecasts (n vector)
  #.........logs........Log-Score of DRN forecasts (n vector)
  #.........lgt.........Length of DRN prediction interval (n vector)
  #.........e_md........Bias of median forecast (n vector)
  #.........e_me........Bias of mean forecast (n vector)
  ###-----------------------------------------------------------------------------
  
  #### Initiation ####
  # Disable GPU
  Sys.setenv("CUDA_VISIBLE_DEVICES" = -1)
  
  # Load packages
  library(keras)
  library(tensorflow)
  
  # Relevant variables for training and testing
  train_vars <- unique(c(obs_var, "location", pred_vars))
  test_vars <- unique(c("location", pred_vars))
  
  # Observations for scores
  if(scores_pp){ test_vars <- c(test_vars, obs_var) }
  
  # Input check
  if(any(!is.element(train_vars, names(train)))){ 
    print("DRN: Training data does not include relevant variables.") }
  if(any(!is.element(test_vars, names(X)))){ 
    print("DRN: Test data does not include all of the relevant variables.") }
  
  # Cut training/test data to relevant variables
  train <- train[,train_vars]
  X <- X[,test_vars]
  
  # Input check
  if(any(is.na(train))){
    print("DRN: Training data includes missing values! Missing values are left out!")
    train <- na.omit(train)
  }
  if(any(is.na(X))){
    print("DRN: Test data includes missing values! Missing values are left out!")
    X <- na.omit(X)
  }
  
  # Number of cores
  if(is.null(n_cores)){ n_cores <- parallel::detectCores()/2 - 1 }
  
  # Set number of cores
  n_cores_config <- tf$compat$v1$ConfigProto(intra_op_parallelism_threads = as.integer(n_cores), 
                                             inter_op_parallelism_threads = as.integer(n_cores))
  tf$compat$v1$keras$backend$set_session(tf$compat$v1$Session(config = n_cores_config))
  
  # Threshold for (almost) constant predictors
  t_c <- 0
  
  #### Hyperparameter ####
  # Hyperparameters and their default values
  hpar_ls <- list(n_sim = 10,
                  distr = "tnorm",
                  t_scale = 1e-4,
                  lr_adam = 5e-4,
                  n_epochs = 150,
                  n_patience = 10,
                  n_batch = 64,
                  emb_dim = 10,
                  layers = c(64, 32),
                  actv = "softplus",
                  nn_verbose = 0)
  
  # Update hyperparameters
  nn_ls <- update_hpar(hpar_ls = hpar_ls,
                       in_ls = nn_ls)
  
  # Custom optimizer
  if(nn_ls$lr_adam == -1){ custom_opt <- "adam" 
  }else{ custom_opt <- optimizer_adam(lr = nn_ls$lr_adam) }
  
  # Get tensorflow probability
  tfp <- tf_probability()$distributions
  
  # Depending on distribution
  if(nn_ls$distr == "tlogis"){
    # Get standard logistic distribution
    tfd <- tfp$Logistic(loc = 0, scale = 1)
    
    # Truncated logistic: Custom loss function
    custom_loss <- function(y_true, y_pred){
      # Get location and scale
      mu  <- k_dot(y_pred, k_constant(as.numeric(rbind(1, 0)), shape = c(2, 1)))
      sigma  <- k_dot(y_pred, k_constant(as.numeric(rbind(0, 1)), shape = c(2, 1)))
      
      # Truncated in zero
      z <- k_maximum(0, y_true)
      
      # Truncate scale to avoid numerical problems
      sigma <- k_maximum(nn_ls$t_scale, sigma)
      
      # Standardization
      z_0 <- -mu/sigma
      z_y <- (z - mu)/sigma
      
      # Calculate CDFs
      p_0 <- tfd$cdf(z_0) # F(z_0)
      lp_0 <- tfd$log_cdf(z_0) # log( F(z_0) )
      p_m0 <- tfd$cdf(-z_0) # 1 - F(z_0) = F(-z_0)
      lp_m0 <- tfd$log_cdf(-z_0) # log( 1 - F(z_0)) = log( F(-z_0) )
      lp_my <- tfd$log_cdf(-z_y) # log( 1 - F(z_y)) = log( F(-z_y) )
      
      # Calculate sigma
      b <- lp_m0 - (1 + 2*lp_my)/p_m0 - k_square(p_0)*lp_0/k_square(p_m0)
      
      # Calculate CRPS
      res <- k_abs(z - y_true) - (z - mu)*(1 + p_0)/p_m0 + sigma*b
      
      # Calculate mean
      res <- k_mean(res)
      
      # Return mean CRPS
      return(res)
    }
  }else if(nn_ls$distr == "tnorm"){
    # Get standard normal distribution
    tfd <- tfp$Normal(loc = 0, scale = 1)
    
    # Truncated normal: Custom loss function (based on formula; Gneiting et al., 2006)
    custom_loss <- function(y_true, y_pred){
      # Get location and scale
      mu  <- k_dot(y_pred, k_constant(as.numeric(rbind(1, 0)), shape = c(2, 1)))
      sigma  <- k_dot(y_pred, k_constant(as.numeric(rbind(0, 1)), shape = c(2, 1)))
      
      # Truncated in zero
      z <- k_maximum(0, y_true)
      
      # Truncate scale to avoid numerical problems
      sigma <- k_maximum(nn_ls$t_scale, sigma)
      
      # Standardization
      z_0 <- mu/sigma
      z_y <- (z - mu)/sigma
      
      # Calculate CDFs and PDF
      p_0 <- tfd$cdf(z_0)
      p_02 <- tfd$cdf(sqrt(2)*z_0)
      p_y <- tfd$cdf(z_y)
      d_y <- tfd$prob(z_y)
  
      # Calculate second term
      b <- z_y*p_0*(2*p_y + p_0 - 2) + 2*d_y*p_0 - p_02/(sqrt(pi))
      
      # Calculate CRPS
      res <- k_abs(z - y_true) + b*sigma/(k_square(p_0))
      
      # Calculate mean
      res <- k_mean(res)
      
      # Return mean CRPS
      return(res)
    }
  }
  
  #### Data preparation ####
  # Divide data in training and validation set
  if(is.null(i_valid)){
    # Set ratio of validation data to 0.2
    r_valid <- 0.2
    
    # Take first n_train*r_valid samples for training
    i_train <- 1:floor(nrow(train)*(1 - r_valid))
    i_valid <- (max(i_train) + 1):nrow(train)
  }
  else{ i_train <- (1:nrow(train))[-i_valid] }
  
  # Read out set sizes
  n_train <- length(i_train)
  n_valid <- length(i_valid)
  n_test <- nrow(X)
  
  # Transform locations to IDs
  train[["location"]] <- sapply(train[["location"]], function(x) which(loc_id_vec == x))
  X[["location"]] <- sapply(X[["location"]], function(x) which(loc_id_vec == x))

  # Remove constant predictors
  pred_vars <- rm_const(data = train[i_train,],
                        cols = pred_vars,
                        t_c = t_c)
  
  # Predictors without location
  dir_preds <- pred_vars[pred_vars != "location"]
  
  # Get number of direct predictors
  n_dir_preds <- length(dir_preds)
  
  # Scale training data of direct predictors
  X_train <- scale(train[i_train, dir_preds])
  
  # Save center and scale parameters
  tr_center <- attr(X_train, "scaled:center")
  tr_scale <- attr(X_train, "scaled:scale")
  
  # Scale validation data with training data attributes
  X_valid <- scale(train[i_valid, dir_preds],
                   center = tr_center,
                   scale = tr_scale)
  
  # Input for fit
  X_train <- list(dir_input = as.matrix(X_train), 
                  id_input = train[i_train, "location"])
  X_valid <- list(as.matrix(X_valid), 
                  train[i_valid, "location"])
  
  # Scale data for prediction
  X_pred <- list(as.matrix(scale(X[,dir_preds], 
                                 center = tr_center, 
                                 scale = tr_scale)), 
                 X[["location"]])
  
  #### Ensemble of networks ####
  # Initiate runtimes
  runtime_est <- runtime_pred <- 0
  
  # Initiate distribution parameter matrix
  f <- matrix(data = 0,
              nrow = nrow(X),
              ncol = 2)
  
  # List for models
  if(model_out){ model_ls <- list() }
  
  # For-Loop over ensemble size
  for(i_sim in 1:nn_ls$n_sim){
    #### Build network ####
    # Input
    id_input <- layer_input(shape = 1, name = "id_input")
    dir_input <- layer_input(shape = n_dir_preds, name = "dir_input")
    
    # Embedding part (see help for input_dim)
    station_embedding_part <- id_input %>%
      layer_embedding(input_dim = length(unique(train[i_train, "location"])) + 1, 
                      output_dim = nn_ls$emb_dim, input_length = 1) %>%
      layer_flatten()
    
    # Concatenate input
    hidden <- layer_concatenate(c(dir_input, station_embedding_part))
    
    # Hidden layers
    for(i_lyr in 1:length(nn_ls$layers)){
      hidden <- hidden %>% 
        layer_dense(units = nn_ls$layers[i_lyr], activation = nn_ls$actv) }
    
    # Different activation functions for output (location positive for numerical reasons)
    loc_out <- layer_dense(object = hidden, units = 1, activation = "softplus")
    scale_out <- layer_dense(object = hidden, units = 1, activation = "softplus")
    
    # Concatenate output
    output <- layer_concatenate(c(loc_out, scale_out))
    
    # Define model
    model <- keras_model(inputs = c(dir_input, id_input), outputs = output)
    
    #### Estimation ####
    # Compile model
    model %>% compile(
      optimizer = custom_opt,
      loss = custom_loss
    )
    
    # Take time
    start_tm <- Sys.time()
    
    # Fit model
    history <- model %>% fit(
      x = X_train,
      y = train[i_train, obs_var],
      epochs = nn_ls$n_epochs,
      batch_size = nn_ls$n_batch,
      validation_data = list(X_valid, train[i_valid, obs_var]),
      verbose = nn_ls$nn_verbose,
      callbacks = callback_early_stopping(patience = nn_ls$n_patience,
                                          restore_best_weights = TRUE,
                                          monitor = "val_loss")
    )

    # Take time
    end_tm <- Sys.time()
    
    # Time needed
    runtime_est <- runtime_est + as.numeric(difftime(end_tm, start_tm, units = "mins"))
    
    # Save model
    if(model_out){ model_ls[[i_sim]] <- model }
    
    # Delete history
    rm(history)
    
    #### Prediction ####
    # Take time
    start_tm <- Sys.time()
    
    # Predict parameters of distributional forecasts (on scaled data)
    f_pred <- predict(model, X_pred)
    
    # Truncate scale parameters
    f_pred[,2] <- pmax(f_pred[,2], nn_ls$t_scale)
    
    # Take time
    end_tm <- Sys.time()
    
    # Time needed
    runtime_pred <- runtime_pred + as.numeric(difftime(end_tm, start_tm, units = "mins"))
    
    # Sum up Predicted parameters for averaging
    f <- f + f_pred
    rm(f_pred)
    
    # Delete model
    rm(id_input, dir_input, station_embedding_part, hidden, 
       loc_out, scale_out, output, model)
    
    # Clear memory and session
    gc()
    k_clear_session()
  }
  
  # Delete data
  rm(X_train, X_valid, X_pred)
  
  # Average distribution parameters
  f <- f/nn_ls$n_sim
  
  #### Evaluation ####
  # Calculate evaluation measures of DRN forecasts
  if(scores_pp){ scores_pp <- fn_scores_distr(f = f,
                                              y = X[[obs_var]],
                                              distr = nn_ls$distr) }
  
  #### Output ####
  # Output list
  res <- list(f = f,
              nn_ls = nn_ls,
              scores_pp = scores_pp, 
              pred_vars = pred_vars,
              n_preds = length(pred_vars),
              n_train = n_train,
              n_valid = n_valid,
              n_test = n_test,
              runtime_est = runtime_est,
              runtime_pred = runtime_pred,
              center = tr_center,
              scale = tr_scale)
  
  # Include model?
  if(model_out){ res[["model_ls"]] <- model_ls }
  
  # Output
  return(res)
}

#### DRN: Bias ####
## Distribution: (Truncated) Normal distribution
## Estimation: CRPS
## Network emsemble aggregation: Parameter averaging (Quantile averaging, Vincentization)
## Comment: Station embedding

# Function for pp including estimation and prediction
drn_pp_bias <- function(train, X, i_valid = NULL, loc_id_vec = NULL,
                        pred_vars = c("ssrd", "location"), 
                        obs_var = "ghi", fc_var = "ssrd", bias_var = "bias",
                        nn_ls = list(), n_cores = NULL, model_out = FALSE,
                        scores_pp = TRUE){
  ###-----------------------------------------------------------------------------
  ###Input
  #train...........Training data including predictors and obs. (n_train x (n_preds + 1) data.frame)
  #X...............Test data including predictors (and obs.) (-"-)
  #i_valid.........Indices of validation data (n_valid vector)
  #................Default: NULL -> Use fraction of training data
  #loc_id_vec......IDs of locations (string vector)
  #................Default: NULL -> Only needed if locations are included
  #pred_vars.......Predictors used for NN (vector of strings)
  #................Default: c("ssrd", "location")
  #obs_var.........Name of observed variable (string)
  #................Default: "ghi"
  #fc_var..........Name of variable to forecast (string)
  #................Default: "ghi"
  #bias_var........Name of bias variable (string)
  #................Default: "bias"
  #nn_ls...........List that may contain the following variables:
  #...n_sim........Number of NN to estimate and average (positive integer)
  #................Default: 10
  #...distr........Forecast distribution ("norm", "logis")
  #................Default: "norm"
  #...t_scale......Threshold for scale parameter truncation (positive scalar)
  #................Default: 1e-4
  #...lr_adam......Learning rate in Adam-Optimizer (scalar)
  #................Default: 5e-4
  #...n_epochs.....Number of epochs (integer)
  #................Default: 150
  #...n_patience...Patience in callbacks (integer)
  #................Default: 10
  #...n_batch......Size of batches is equal to n_train/n_batch (input parameter batch_size)
  #................Default: 32
  #...emb_dim......Dimension of station embedding (scalar)
  #................Default: 5
  #...layers.......Vector that defines the nodes in the hidden layers (n_layers vector)
  #................Default: c(96, 48, 24) -> 3 layers
  #...actv.........Activation function of non-output layers (string)
  #................Default: "softplus"
  #...nn_verbose...Query, if network output should be shown (logical)
  #................Default: 0
  #n_cores.........Number of cores used in keras (integer)
  #................Default: NULL -> Use one less than available
  #model_out.......Query, if models should be in output (logical)
  #................Default: FALSE -> No model output
  #scores_pp.......Should evaluation measures be calculated (logical)
  #................Default: TRUE
  ###-----------------------------------------------------------------------------
  ###Output
  #res...List containing:
  #......f..............Distributional forecasts (n x 2 matrix)
  #......nn_ls..........Hyperparameters (list)
  #......model_ls.......Network models of the ensemble (list)
  #......pred_vars......Predictors (string vector)
  #......n_preds........Number of predictors (integer)
  #......n_train........Number of training samples (integer)
  #......n_valid........Number of validation samples (integer)
  #......n_test.........Number of test samples (integer)
  #......runtime_est....Estimation time (numeric)
  #......runtime_pred...Prediction time (numeric)
  #......center.........Centering parameter of training data (numeric)
  #......scale..........Scaling parameter of training data (numeric)
  #......scores_pp......Data frames containing (n x 6 data frame):
  #.........pit.........PIT of DRN forecasts (n vector)
  #.........crps........CRPS of DRN forecasts (n vector)
  #.........logs........Log-Score of DRN forecasts (n vector)
  #.........lgt.........Length of DRN prediction interval (n vector)
  #.........e_md........Bias of median forecast (n vector)
  #.........e_me........Bias of mean forecast (n vector)
  ###-----------------------------------------------------------------------------
  
  #### Initiation ####
  # Disable GPU
  Sys.setenv("CUDA_VISIBLE_DEVICES" = -1)
  
  # Load packages
  library(keras)
  library(tensorflow)
  
  # Relevant variables for training and testing
  train_vars <- unique(c(bias_var, "location", pred_vars))
  test_vars <- unique(c("location", pred_vars, fc_var, bias_var))
  
  # Observations for scores
  if(scores_pp){ test_vars <- c(test_vars, obs_var) }
  
  # Input check
  if(any(!is.element(train_vars, names(train)))){ 
    print("DRN: Training data does not include relevant variables.") }
  if(any(!is.element(test_vars, names(X)))){ 
    print("DRN: Test data does not include all of the relevant variables.") }
  
  # Cut training/test data to relevant variables
  train <- train[,train_vars]
  X <- X[,test_vars]
  
  # Input check
  if(any(is.na(train))){
    print("DRN: Training data includes missing values! Missing values are left out!")
    train <- na.omit(train)
  }
  if(any(is.na(X))){
    print("DRN: Test data includes missing values! Missing values are left out!")
    X <- na.omit(X)
  }
  
  # Number of cores
  if(is.null(n_cores)){ n_cores <- parallel::detectCores()/2 - 1 }
  
  # Set number of cores
  n_cores_config <- tf$compat$v1$ConfigProto(intra_op_parallelism_threads = as.integer(n_cores), 
                                             inter_op_parallelism_threads = as.integer(n_cores))
  tf$compat$v1$keras$backend$set_session(tf$compat$v1$Session(config = n_cores_config))
  
  # Threshold for (almost) constant predictors
  t_c <- 0
  
  #### Hyperparameter ####
  # Hyperparameters and their default values
  hpar_ls <- list(n_sim = 10,
                  distr = "norm",
                  t_scale = 1e-4,
                  lr_adam = 5e-4,
                  n_epochs = 150,
                  n_patience = 10,
                  n_batch = 32,
                  emb_dim = 5,
                  layers = c(96, 48, 24),
                  actv = "softplus",
                  nn_verbose = 0)
  
  # Update hyperparameters
  nn_ls <- update_hpar(hpar_ls = hpar_ls,
                       in_ls = nn_ls)
  
  # Custom optimizer
  if(nn_ls$lr_adam == -1){ custom_opt <- "adam" 
  }else{ custom_opt <- optimizer_adam(lr = nn_ls$lr_adam) }
  
  # Get tensorflow probability
  tfp <- tf_probability()$distributions
  
  # Depending on distribution
  if(nn_ls$distr == "norm"){
    # Get standard normal distribution
    tfd <- tfp$Normal(loc = 0, scale = 1)
    
    # Normal: Custom loss function
    custom_loss <- function(y_true, y_pred){
      # Get location and scale
      mu  <- k_dot(y_pred, k_constant(as.numeric(rbind(1, 0)), shape = c(2, 1)))
      sigma  <- k_dot(y_pred, k_constant(as.numeric(rbind(0, 1)), shape = c(2, 1)))
      
      # Truncate scale to avoid numerical problems
      sigma <- k_maximum(nn_ls$t_scale, sigma)
      
      # Standardization
      z_y <- (y_true - mu)/sigma
      
      # Calculate CRPS
      res <- sigma*( z_y*(2*tfd$cdf(z_y) - 1) + 2*tfd$prob(z_y) - 1/sqrt(pi) )
      
      # Calculate mean
      res <- k_mean(res)
      
      # Return mean CRPS
      return(res)
    }
  }else if(nn_ls$distr == "logis"){
    # Get standard logistic distribution
    tfd <- tfp$Logistic(loc = 0, scale = 1)
    
    # Logistic: Custom loss function
    custom_loss <- function(y_true, y_pred){
      # Get location and scale
      mu  <- k_dot(y_pred, k_constant(as.numeric(rbind(1, 0)), shape = c(2, 1)))
      sigma  <- k_dot(y_pred, k_constant(as.numeric(rbind(0, 1)), shape = c(2, 1)))
      
      # Truncate scale to avoid numerical problems
      sigma <- k_maximum(nn_ls$t_scale, sigma)
      
      # Standardization
      z_y <- (y_true - mu)/sigma
      
      # Calculate CRPS
      res <- sigma*( z_y - 2*tfd$log_cdf(z_y) - 1 )

      # Calculate mean
      res <- k_mean(res)
      
      # Return mean CRPS
      return(res)
    }
  }
  
  #### Data preparation ####
  # Divide data in training and validation set
  if(is.null(i_valid)){
    # Set ratio of validation data to 0.2
    r_valid <- 0.2
    
    # Take first n_train*r_valid samples for training
    i_train <- 1:floor(nrow(train)*(1 - r_valid))
    i_valid <- (max(i_train) + 1):nrow(train)
  }
  else{ i_train <- (1:nrow(train))[-i_valid] }
  
  # Read out set sizes
  n_train <- length(i_train)
  n_valid <- length(i_valid)
  n_test <- nrow(X)
  
  # Transform locations to IDs
  train[["location"]] <- sapply(train[["location"]], function(x) which(loc_id_vec == x))
  X[["location"]] <- sapply(X[["location"]], function(x) which(loc_id_vec == x))
  
  # Remove constant predictors
  pred_vars <- rm_const(data = train[i_train,],
                        cols = pred_vars,
                        t_c = t_c)
  
  # Predictors without location
  dir_preds <- pred_vars[pred_vars != "location"]
  
  # Get number of direct predictors
  n_dir_preds <- length(dir_preds)
  
  # Scale training data of direct predictors
  X_train <- scale(train[i_train, dir_preds])
  
  # Save center and scale parameters
  tr_center <- attr(X_train, "scaled:center")
  tr_scale <- attr(X_train, "scaled:scale")
  
  # Scale validation data with training data attributes
  X_valid <- scale(train[i_valid, dir_preds],
                   center = tr_center,
                   scale = tr_scale)
  
  # Input for fit
  X_train <- list(dir_input = as.matrix(X_train), 
                  id_input = train[i_train, "location"])
  X_valid <- list(as.matrix(X_valid), 
                  train[i_valid, "location"])
  
  # Scale data for prediction
  X_pred <- list(as.matrix(scale(X[,dir_preds], 
                                 center = tr_center, 
                                 scale = tr_scale)), 
                 X[["location"]])
  
  #### Ensemble of networks ####
  # Initiate runtimes
  runtime_est <- runtime_pred <- 0
  
  # Initiate distribution parameter matrix
  f <- matrix(data = 0,
              nrow = nrow(X),
              ncol = 2)
  
  # List for models
  if(model_out){ model_ls <- list() }
  
  # For-Loop over ensemble size
  for(i_sim in 1:nn_ls$n_sim){
    #### Build network ####
    # Input
    id_input <- layer_input(shape = 1, name = "id_input")
    dir_input <- layer_input(shape = n_dir_preds, name = "dir_input")
    
    # Embedding part (see help for input_dim)
    station_embedding_part <- id_input %>%
      layer_embedding(input_dim = length(unique(train[i_train, "location"])) + 1, 
                      output_dim = nn_ls$emb_dim, input_length = 1) %>%
      layer_flatten()
    
    # Concatenate input
    hidden <- layer_concatenate(c(dir_input, station_embedding_part))
    
    # Hidden layers
    for(i_lyr in 1:length(nn_ls$layers)){
      hidden <- hidden %>% 
        layer_dense(units = nn_ls$layers[i_lyr], activation = nn_ls$actv) }
    
    # Different activation functions for output
    loc_out <- layer_dense(object = hidden, units = 1)
    scale_out <- layer_dense(object = hidden, units = 1, activation = "softplus")
    
    # Concatenate output
    output <- layer_concatenate(c(loc_out, scale_out))
    
    # Define model
    model <- keras_model(inputs = c(dir_input, id_input), outputs = output)
    
    #### Estimation ####
    # Compile model
    model %>% compile(
      optimizer = custom_opt,
      loss = custom_loss
    )
    
    # Take time
    start_tm <- Sys.time()
    
    # Fit model
    history <- model %>% fit(
      x = X_train,
      y = train[i_train, bias_var],
      epochs = nn_ls$n_epochs,
      batch_size = nn_ls$n_batch,
      validation_data = list(X_valid, train[i_valid, bias_var]),
      verbose = nn_ls$nn_verbose,
      callbacks = callback_early_stopping(patience = nn_ls$n_patience,
                                          restore_best_weights = TRUE,
                                          monitor = "val_loss")
    )
    
    # Take time
    end_tm <- Sys.time()
    
    # Time needed
    runtime_est <- runtime_est + as.numeric(difftime(end_tm, start_tm, units = "mins"))
    
    # Save model
    if(model_out){ model_ls[[i_sim]] <- model }
    
    # Delete history
    rm(history)
    
    #### Prediction ####
    # Take time
    start_tm <- Sys.time()
    
    # Predict parameters of distributional forecasts (on scaled data)
    f_pred <- predict(model, X_pred)
    
    # Truncate scale parameters
    f_pred[,2] <- pmax(f_pred[,2], nn_ls$t_scale)
    
    # Take time
    end_tm <- Sys.time()
    
    # Time needed
    runtime_pred <- runtime_pred + as.numeric(difftime(end_tm, start_tm, units = "mins"))
    
    # Sum up Predicted parameters for averaging
    f <- f + f_pred
    rm(f_pred)
    
    # Delete model
    rm(id_input, dir_input, station_embedding_part, hidden, 
       loc_out, scale_out, output, model)
    
    # Clear memory and session
    gc()
    k_clear_session()
  }
  
  # Delete data
  rm(X_train, X_valid, X_pred)
  
  # Average distribution parameters
  f <- f/nn_ls$n_sim
  
  # Add up forecast and bias
  f[,1] <- X[,fc_var] - f[,1] 
  
  #### Evaluation ####
  # Calculate evaluation measures of DRN forecasts 
  # Use a (zero-)truncated distribution for forecasting (!)
  if(scores_pp){ scores_pp <- fn_scores_distr(f = f,
                                              y = X[[obs_var]],
                                              distr = paste0("t", nn_ls$distr)) }
  
  #### Output ####
  # Output list
  res <- list(f = f,
              nn_ls = nn_ls,
              scores_pp = scores_pp, 
              pred_vars = pred_vars,
              n_preds = length(pred_vars),
              n_train = n_train,
              n_valid = n_valid,
              n_test = n_test,
              runtime_est = runtime_est,
              runtime_pred = runtime_pred,
              center = tr_center,
              scale = tr_scale)
  
  # Include model?
  if(model_out){ res[["model_ls"]] <- model_ls }
  
  # Output
  return(res)
}

#### BQN ####
## Network emsemble aggregation: Parameter averaging (Quantile averaging, Vincentization)
## Comment: Station embedding

# Function for pp including estimation and prediction
bqn_pp <- function(train, X, i_valid = NULL, loc_id_vec = NULL,
                   pred_vars = c("ssrd", "location"), obs_var = "ghi",
                   nn_ls = list(), n_cores = NULL, model_out = FALSE,
                   q_levels = NULL, scores_pp = TRUE){
  ###-----------------------------------------------------------------------------
  ###Input
  #train...........Training data including predictors and obs. (n_train x (n_preds + 1) data.frame)
  #X...............Test data including predictors (and obs.) (-"-)
  #i_valid.........Indices of validation data (n_valid vector)
  #................Default: NULL -> Use fraction of training data
  #loc_id_vec......IDs of locations (string vector)
  #................Default: NULL -> Only needed if locations are included
  #pred_vars.......Predictors used for NN (vector of strings)
  #................Default: c("ssrd", "location")
  #obs_var.........Name of observed variable (string)
  #................Default: "ghi"
  #nn_ls...........List that may contain the following variables:
  #...p_degree.....Degree of Bernstein polynomials (integer)
  #................Default: 12
  #...n_q..........Number of equidistant quantile levels used in loss function (integer)
  #................Default: 99 (steps of 1%)
  #...n_sim........Number of NN to estimate and average (positive integer)
  #................Default: 10
  #...lr_adam......Learning rate in Adam-Optimizer (scalar)
  #................Default: 5e-4
  #...n_epochs.....Number of epochs (integer)
  #................Default: 150
  #...n_patience...Patience in callbacks (integer)
  #................Default: 10
  #...n_batch......Size of batches is equal to n_train/n_batch (input parameter batch_size)
  #................Default: 64
  #...emb_dim......Dimension of station embedding (scalar)
  #................Default: 10
  #...layers.......Vector that defines the nodes in the hidden layers (n_layers vector)
  #................Default: c(48, 24) -> 2 layers: 48 in first, 24 in second
  #...actv.........Activation function of non-output layers (string)
  #................Default: "softplus"
  #...actv_out.....Activation function of output layer (string)
  #................Default: "softplus"
  #...nn_verbose...Query, if network output should be shown (logical)
  #................Default: 0
  #n_cores.........Number of cores used in keras (integer)
  #................Default: NULL -> Use one less than available
  #model_out.......Query, if models should be in output (logical)
  #................Default: FALSE -> No model output
  #q_levels........Quantile levels used for output and evaluation (probability vector)
  #................Default: NULL -> 100 member
  #scores_pp.......Should evaluation measures be calculated (logical)
  #................Default: TRUE
  ###-----------------------------------------------------------------------------
  ###Output
  #res...List containing:
  #......f..............BQN quantile forecasts (n x n_q matrix)
  #......q_levels.......Quantile levels of forecasts in f (n_q probability vector)
  #......alpha..........BQN coefficients (n x p_degree matrix)
  #......nn_ls..........Hyperparameters (list)
  #......model_ls.......Network models of the ensemble (list)
  #......pred_vars......Predictors (string vector)
  #......n_preds........Number of predictors (integer)
  #......n_train........Number of training samples (integer)
  #......n_valid........Number of validation samples (integer)
  #......n_test.........Number of test samples (integer)
  #......runtime_est....Estimation time (numeric)
  #......runtime_pred...Prediction time (numeric)
  #......center.........Centering parameter of training data (numeric)
  #......scale..........Scaling parameter of training data (numeric)
  #......scores_pp......Data frames containing (n x 6 data frame):
  #.........pit.........uPIT of BQN forecasts (n vector)
  #.........crps........CRPS of BQN forecasts (n vector)
  #.........logs........Log-Score of BQN forecasts (n vector)
  #.........lgt.........Length of BQN prediction interval (n vector)
  #.........e_md........Bias of median forecast (n vector)
  #.........e_me........Bias of mean forecast (n vector)
  ###-----------------------------------------------------------------------------
  
  #### Initiation ####
  # Disable GPU
  Sys.setenv("CUDA_VISIBLE_DEVICES" = -1)
  
  # Load packages
  library(keras)
  library(tensorflow)
  
  # Relevant variables for training and testing
  train_vars <- unique(c(obs_var, "location", pred_vars))
  test_vars <- unique(c("location", pred_vars))
  
  # Observations for scores
  if(scores_pp){ test_vars <- c(test_vars, obs_var) }
  
  # Input check
  if(any(!is.element(train_vars, names(train)))){ 
    print("BQN: Training data does not include relevant variables.") }
  if(any(!is.element(test_vars, names(X)))){ 
    print("BQN: Test data does not include all of the relevant variables.") }
  
  # Cut training/test data to relevant variables
  train <- train[,train_vars]
  X <- X[,test_vars]
  
  # Input check
  if(any(is.na(train))){
    print("BQN: Training data includes missing values! Missing values are left out!")
    train <- na.omit(train)
  }
  if(any(is.na(X))){
    print("BQN: Test data includes missing values! Missing values are left out!")
    X <- na.omit(X)
  }
  
  # Number of cores
  if(is.null(n_cores)){ n_cores <- parallel::detectCores()/2 - 1 }
  
  # Set number of cores
  n_cores_config <- tf$compat$v1$ConfigProto(intra_op_parallelism_threads = as.integer(n_cores), 
                                             inter_op_parallelism_threads = as.integer(n_cores))
  tf$compat$v1$keras$backend$set_session(tf$compat$v1$Session(config = n_cores_config))
  
  # Threshold for (almost) constant predictors
  t_c <- 0
  
  # If not given use equidistant quantiles (multiple of ensemble coverage, incl. median)
  if(is.null(q_levels)){ q_levels <- seq(from = 1/100,
                                         to = 99/100,
                                         by = 1/100) }
  
  #### Hyperparameter ####
  # Hyperparameters and their default values
  hpar_ls <- list(p_degree = 12,
                  n_q = 99,
                  n_sim = 10,
                  lr_adam = 5e-4, # -1 for Adam-default
                  n_epochs = 150,
                  n_patience = 10,
                  n_batch = 64,
                  emb_dim = 10,
                  layers = c(48, 24),
                  actv = "softplus",
                  actv_out = "softplus",
                  nn_verbose = 0)
  
  # Update hyperparameters
  nn_ls <- update_hpar(hpar_ls = hpar_ls,
                       in_ls = nn_ls)
  
  # Calculate equidistant quantile levels for loss function
  q_levels_loss <- seq(from = 1/(nn_ls$n_q + 1),
                       to = nn_ls$n_q/(nn_ls$n_q + 1),
                       by = 1/(nn_ls$n_q + 1))
  
  # Basis of Bernstein polynomials evaluated at quantile levels
  B <- sapply(0:nn_ls$p_degree, 
              dbinom, size = nn_ls$p_degree, prob = q_levels_loss)
  
  # Quantile loss functions (for neural network)
  qt_loss <- function(y_true, y_pred){
    # Quantiles calculated via basis and increments
    q  <- k_dot(k_cumsum(y_pred, axis = 0),
                k_constant(as.numeric(B), shape = c(nn_ls$p_degree + 1, nn_ls$n_q)))
    
    # Calculate individual quantile scores
    err  <- y_true - q
    e1   <- err * k_constant(q_levels_loss, shape = c(1, nn_ls$n_q))
    e2   <- err * k_constant(q_levels_loss - 1, shape = c(1, nn_ls$n_q))
    
    # Find correct values (max) and return mean
    return(k_mean( k_maximum(e1, e2), axis = 2 ))
  }
  
  # Custom optimizer
  if(nn_ls$lr_adam == -1){ custom_opt <- "adam" 
  }else{ custom_opt <- optimizer_adam(lr = nn_ls$lr_adam) }
  
  #### Data preparation ####
  # Divide data in training and validation set
  if(is.null(i_valid)){
    # Set ratio of validation data to 0.2
    r_valid <- 0.2
    
    # Take first n_train*r_valid samples for training
    i_train <- 1:floor(nrow(train)*(1 - r_valid))
    i_valid <- (max(i_train) + 1):nrow(train)
  }
  else{ i_train <- (1:nrow(train))[-i_valid] }
  
  # Read out set sizes
  n_train <- length(i_train)
  n_valid <- length(i_valid)
  n_test <- nrow(X)
  
  # Transform locations to IDs
  train[["location"]] <- sapply(train[["location"]], function(x) which(loc_id_vec == x))
  X[["location"]] <- sapply(X[["location"]], function(x) which(loc_id_vec == x))
  
  # Remove constant predictors
  pred_vars <- rm_const(data = train[i_train,],
                        cols = pred_vars,
                        t_c = t_c)
  
  # Predictors without location
  dir_preds <- pred_vars[pred_vars != "location"]
  
  # Get number of direct predictors
  n_dir_preds <- length(dir_preds)
  
  # Scale training data of direct predictors
  X_train <- scale(train[i_train, dir_preds])
  
  # Save center and scale parameters
  tr_center <- attr(X_train, "scaled:center")
  tr_scale <- attr(X_train, "scaled:scale")
  
  # Scale validation data with training data attributes
  X_valid <- scale(train[i_valid, dir_preds],
                   center = tr_center,
                   scale = tr_scale)
  
  # Input for fit
  X_train <- list(dir_input = as.matrix(X_train), 
                  id_input = train[i_train, "location"])
  X_valid <- list(as.matrix(X_valid), 
                  train[i_valid, "location"])

  # Scale data for prediction
  X_pred <- list(as.matrix(scale(X[,dir_preds], 
                                 center = tr_center, 
                                 scale = tr_scale)), 
                 X[["location"]])
  
  #### Ensemble of networks ####
  # Initiate runtimes
  runtime_est <- runtime_pred <- 0
  
  # Initiate coefficient matrix
  coeff_bern <- matrix(data = 0,
                       nrow = nrow(X),
                       ncol = nn_ls$p_degree + 1)
  
  # List for models
  if(model_out){ model_ls <- list() }
  
  # For-Loop over ensemble size
  for(i_sim in 1:nn_ls$n_sim){
    #### Build network ####
    # Input
    id_input <- layer_input(shape = 1, name = "id_input")
    dir_input <- layer_input(shape = n_dir_preds, name = "dir_input")
    
    # Embedding part (see help for input_dim)
    station_embedding_part <- id_input %>%
      layer_embedding(input_dim = length(unique(train[i_train, "location"])) + 1, 
                      output_dim = nn_ls$emb_dim, input_length = 1) %>%
      layer_flatten()
    
    # Concatenate input
    hidden <- layer_concatenate(c(dir_input, station_embedding_part))
    
    # Hidden layers
    for(i_lyr in 1:length(nn_ls$layers)){
      hidden <- hidden %>% 
        layer_dense(units = nn_ls$layers[i_lyr], activation = nn_ls$actv) }
    
    # Monotonicity in output (increments; strict via softplus, else relu)
    output <- hidden %>%
      layer_dense(units = nn_ls$p_degree + 1, activation = nn_ls$actv_out)
    
    # Define model
    model <- keras_model(inputs = c(dir_input, id_input), outputs = output)
    
    #### Estimation ####
    # Compile model
    model %>% compile(
      optimizer = custom_opt,
      loss = qt_loss
    )
    
    # Take time
    start_tm <- Sys.time()
    
    # Fit model
    history <- model %>% fit(
      x = X_train,
      y = train[i_train, obs_var],
      epochs = nn_ls$n_epochs,
      batch_size = nn_ls$n_batch,
      validation_data = list(X_valid, train[i_valid, obs_var]),
      verbose = nn_ls$nn_verbose,
      callbacks = callback_early_stopping(patience = nn_ls$n_patience,
                                          restore_best_weights = TRUE,
                                          monitor = "val_loss")
    )
    
    # Take time
    end_tm <- Sys.time()
    
    # Time needed
    runtime_est <- runtime_est + as.numeric(difftime(end_tm, start_tm, units = "mins"))
    
    # Save model
    if(model_out){ model_ls[[i_sim]] <- model }
    
    # Delete history
    rm(history)
    
    #### Prediction ####
    # Take time
    start_tm <- Sys.time()
    
    # Predict coefficients of Bernstein polynomials
    coeff_bern <- coeff_bern + predict(model, X_pred)
    
    # Take time
    end_tm <- Sys.time()
    
    # Time needed
    runtime_pred <- runtime_pred + as.numeric(difftime(end_tm, start_tm, units = "mins"))
    
    # Delete model
    rm(id_input, dir_input, station_embedding_part, hidden, 
       output, model)
    
    # Clear memory and session
    gc()
    k_clear_session()
  }
  
  # Delete data
  rm(X_train, X_valid, X_pred)
  
  # Average increments
  coeff_bern <- coeff_bern/nn_ls$n_sim
  
  # Accumulate increments
  coeff_bern <- t(apply(coeff_bern, 1, cumsum))
  
  #### Evaluation ####
  # Sum up calculated quantiles (Sum of basis at quantiles times coefficients)
  q <- bern_quants(alpha = coeff_bern,
                   q_levels = q_levels)
  
  # Calculate evaluation measure of BQN forecasts
  scores_pp <- fn_scores_ens(ens = q,
                             y = X[[obs_var]],
                             skip_evals = c("e_me"),
                             scores_ens = scores_pp)
  
  # Transform ranks to uPIT
  scores_pp[["pit"]] <- scores_pp[["rank"]]/(ncol(q) + 1) - 
    runif(n = n_test,
          min = 0,
          max = 1/(ncol(q) + 1))
  scores_pp[["rank"]] <- NULL
  
  # Calculate bias of mean forecast (formula given)
  scores_pp[["e_me"]] <- rowMeans(coeff_bern) - X[[obs_var]]
  
  #### Output ####
  # Output list
  res <- list(f = q, 
              q_levels = q_levels,
              alpha = coeff_bern,
              nn_ls = nn_ls,
              scores_pp = scores_pp, 
              pred_vars = pred_vars,
              n_preds = length(pred_vars),
              n_train = n_train,
              n_valid = n_valid,
              n_test = n_test,
              runtime_est = runtime_est,
              runtime_pred = runtime_pred,
              center = tr_center,
              scale = tr_scale)
  
  # Include model?
  if(model_out){ res[["model_ls"]] <- model_ls }
  
  # Output
  return(res)
}

#### BQN: Bias ####
## Network emsemble aggregation: Parameter averaging (Quantile averaging, Vincentization)
## Comment: Station embedding

# Function for pp including estimation and prediction
bqn_pp_bias <- function(train, X, i_valid = NULL, loc_id_vec = NULL,
                   pred_vars = c("ssrd", "location"), 
                   obs_var = "ghi", fc_var = "ssrd", bias_var = "bias",
                   nn_ls = list(), n_cores = NULL, model_out = FALSE,
                   q_levels = NULL, scores_pp = TRUE){
  ###-----------------------------------------------------------------------------
  ###Input
  #train...........Training data including predictors and obs. (n_train x (n_preds + 1) data.frame)
  #X...............Test data including predictors (and obs.) (-"-)
  #i_valid.........Indices of validation data (n_valid vector)
  #................Default: NULL -> Use fraction of training data
  #loc_id_vec......IDs of locations (string vector)
  #................Default: NULL -> Only needed if locations are included
  #pred_vars.......Predictors used for NN (vector of strings)
  #................Default: c("ssrd", "location")
  #obs_var.........Name of observed variable (string)
  #................Default: "ghi"
  #fc_var..........Name of variable to forecast (string)
  #................Default: "ghi"
  #bias_var........Name of bias variable (string)
  #................Default: "bias"
  #nn_ls...........List that may contain the following variables:
  #...p_degree.....Degree of Bernstein polynomials (integer)
  #................Default: 12
  #...n_q..........Number of equidistant quantile levels used in loss function (integer)
  #................Default: 99 (steps of 1%)
  #...n_sim........Number of NN to estimate and average (positive integer)
  #................Default: 10
  #...lr_adam......Learning rate in Adam-Optimizer (scalar)
  #................Default: 5e-4
  #...n_epochs.....Number of epochs (integer)
  #................Default: 150
  #...n_patience...Patience in callbacks (integer)
  #................Default: 10
  #...n_batch......Size of batches is equal to n_train/n_batch (input parameter batch_size)
  #................Default: 32
  #...emb_dim......Dimension of station embedding (scalar)
  #................Default: 5
  #...layers.......Vector that defines the nodes in the hidden layers (n_layers vector)
  #................Default: c(96, 48, 24) -> 3 layers
  #...actv.........Activation function of non-output layers (string)
  #................Default: "softplus"
  #...actv_out.....Activation function of output layer (string)
  #................Default: "softplus"
  #...nn_verbose...Query, if network output should be shown (logical)
  #................Default: 0
  #n_cores.........Number of cores used in keras (integer)
  #................Default: NULL -> Use one less than available
  #model_out.......Query, if models should be in output (logical)
  #................Default: FALSE -> No model output
  #q_levels........Quantile levels used for output and evaluation (probability vector)
  #................Default: NULL -> 99 member from 0.01 to 0.99
  #scores_pp.......Should evaluation measures be calculated (logical)
  #................Default: TRUE
  ###-----------------------------------------------------------------------------
  ###Output
  #res...List containing:
  #......f..............BQN quantile forecasts (n x n_q matrix)
  #......q_levels.......Quantile levels of forecasts in f (n_q probability vector)
  #......alpha..........BQN coefficients (n x p_degree matrix)
  #......nn_ls..........Hyperparameters (list)
  #......model_ls.......Network models of the ensemble (list)
  #......pred_vars......Predictors (string vector)
  #......n_preds........Number of predictors (integer)
  #......n_train........Number of training samples (integer)
  #......n_valid........Number of validation samples (integer)
  #......n_test.........Number of test samples (integer)
  #......runtime_est....Estimation time (numeric)
  #......runtime_pred...Prediction time (numeric)
  #......center.........Centering parameter of training data (numeric)
  #......scale..........Scaling parameter of training data (numeric)
  #......scores_pp......Data frames containing (n x 6 data frame):
  #.........pit.........uPIT of BQN forecasts (n vector)
  #.........crps........CRPS of BQN forecasts (n vector)
  #.........logs........Log-Score of BQN forecasts (n vector)
  #.........lgt.........Length of BQN prediction interval (n vector)
  #.........e_md........Bias of median forecast (n vector)
  #.........e_me........Bias of mean forecast (n vector)
  ###-----------------------------------------------------------------------------
  
  #### Initiation ####
  # Disable GPU
  Sys.setenv("CUDA_VISIBLE_DEVICES" = -1)
  
  # Load packages
  library(keras)
  library(tensorflow)
  
  # Relevant variables for training and testing
  train_vars <- unique(c(bias_var, "location", pred_vars))
  test_vars <- unique(c("location", pred_vars, fc_var, bias_var))
  
  # Observations for scores
  if(scores_pp){ test_vars <- c(test_vars, obs_var) }
  
  # Input check
  if(any(!is.element(train_vars, names(train)))){ 
    print("BQN: Training data does not include relevant variables.") }
  if(any(!is.element(test_vars, names(X)))){ 
    print("BQN: Test data does not include all of the relevant variables.") }
  
  # Cut training/test data to relevant variables
  train <- train[,train_vars]
  X <- X[,test_vars]
  
  # Input check
  if(any(is.na(train))){
    print("BQN: Training data includes missing values! Missing values are left out!")
    train <- na.omit(train)
  }
  if(any(is.na(X))){
    print("BQN: Test data includes missing values! Missing values are left out!")
    X <- na.omit(X)
  }
  
  # Number of cores
  if(is.null(n_cores)){ n_cores <- parallel::detectCores()/2 - 1 }
  
  # Set number of cores
  n_cores_config <- tf$compat$v1$ConfigProto(intra_op_parallelism_threads = as.integer(n_cores), 
                                             inter_op_parallelism_threads = as.integer(n_cores))
  tf$compat$v1$keras$backend$set_session(tf$compat$v1$Session(config = n_cores_config))
  
  # Threshold for (almost) constant predictors
  t_c <- 0
  
  # If not given use equidistant quantiles (multiple of ensemble coverage, incl. median)
  if(is.null(q_levels)){ q_levels <- seq(from = 1/100,
                                         to = 99/100,
                                         by = 1/100) }
  
  #### Hyperparameter ####
  # Hyperparameters and their default values
  hpar_ls <- list(p_degree = 12,
                  n_q = 99, 
                  n_sim = 10,
                  lr_adam = 5e-4, # -1 for Adam-default
                  n_epochs = 150,
                  n_patience = 10,
                  n_batch = 32,
                  emb_dim = 5,
                  layers = c(96, 48, 24),
                  actv = "softplus",
                  actv_out = "softplus",
                  nn_verbose = 0)
  
  # Update hyperparameters
  nn_ls <- update_hpar(hpar_ls = hpar_ls,
                       in_ls = nn_ls)
  
  # Calculate equidistant quantile levels for loss function
  q_levels_loss <- seq(from = 1/(nn_ls$n_q + 1),
                       to = nn_ls$n_q/(nn_ls$n_q + 1),
                       by = 1/(nn_ls$n_q + 1))
  
  # Basis of Bernstein polynomials evaluated at quantile levels
  B <- sapply(0:nn_ls$p_degree, 
              dbinom, size = nn_ls$p_degree, prob = q_levels_loss)
  
  # Quantile loss functions (for neural network)
  qt_loss <- function(y_true, y_pred){
    # Quantiles calculated via basis and increments
    q  <- k_dot(k_cumsum(y_pred, axis = 0),
                k_constant(as.numeric(B), shape = c(nn_ls$p_degree + 1, nn_ls$n_q)))
    
    # Calculate individual quantile scores
    err  <- y_true - q
    e1   <- err * k_constant(q_levels_loss, shape = c(1, nn_ls$n_q))
    e2   <- err * k_constant(q_levels_loss - 1, shape = c(1, nn_ls$n_q))
    
    # Find correct values (max) and return mean
    return(k_mean( k_maximum(e1, e2), axis = 2 ))
  }
  
  # Custom optimizer
  if(nn_ls$lr_adam == -1){ custom_opt <- "adam" 
  }else{ custom_opt <- optimizer_adam(lr = nn_ls$lr_adam) }
  
  #### Data preparation ####
  # Divide data in training and validation set
  if(is.null(i_valid)){
    # Set ratio of validation data to 0.2
    r_valid <- 0.2
    
    # Take first n_train*r_valid samples for training
    i_train <- 1:floor(nrow(train)*(1 - r_valid))
    i_valid <- (max(i_train) + 1):nrow(train)
  }
  else{ i_train <- (1:nrow(train))[-i_valid] }
  
  # Read out set sizes
  n_train <- length(i_train)
  n_valid <- length(i_valid)
  n_test <- nrow(X)
  
  # Transform locations to IDs
  train[["location"]] <- sapply(train[["location"]], function(x) which(loc_id_vec == x))
  X[["location"]] <- sapply(X[["location"]], function(x) which(loc_id_vec == x))
  
  # Remove constant predictors
  pred_vars <- rm_const(data = train[i_train,],
                        cols = pred_vars,
                        t_c = t_c)
  
  # Predictors without location
  dir_preds <- pred_vars[pred_vars != "location"]
  
  # Get number of direct predictors
  n_dir_preds <- length(dir_preds)
  
  # Scale training data of direct predictors
  X_train <- scale(train[i_train, dir_preds])
  
  # Save center and scale parameters
  tr_center <- attr(X_train, "scaled:center")
  tr_scale <- attr(X_train, "scaled:scale")
  
  # Scale validation data with training data attributes
  X_valid <- scale(train[i_valid, dir_preds],
                   center = tr_center,
                   scale = tr_scale)
  
  # Input for fit
  X_train <- list(dir_input = as.matrix(X_train), 
                  id_input = train[i_train, "location"])
  X_valid <- list(as.matrix(X_valid), 
                  train[i_valid, "location"])
  
  # Scale data for prediction
  X_pred <- list(as.matrix(scale(X[,dir_preds], 
                                 center = tr_center, 
                                 scale = tr_scale)), 
                 X[["location"]])
  
  #### Ensemble of networks ####
  # Initiate runtimes
  runtime_est <- runtime_pred <- 0
  
  # Initiate coefficient matrix
  coeff_bern <- matrix(data = 0,
                       nrow = nrow(X),
                       ncol = nn_ls$p_degree + 1)
  
  # List for models
  if(model_out){ model_ls <- list() }
  
  # For-Loop over ensemble size
  for(i_sim in 1:nn_ls$n_sim){
    #### Build network ####
    # Input
    id_input <- layer_input(shape = 1, name = "id_input")
    dir_input <- layer_input(shape = n_dir_preds, name = "dir_input")
    
    # Embedding part (see help for input_dim)
    station_embedding_part <- id_input %>%
      layer_embedding(input_dim = length(unique(train[i_train, "location"])) + 1, 
                      output_dim = nn_ls$emb_dim, input_length = 1) %>%
      layer_flatten()
    
    # Concatenate input
    hidden <- layer_concatenate(c(dir_input, station_embedding_part))
    
    # Hidden layers
    for(i_lyr in 1:length(nn_ls$layers)){
      hidden <- hidden %>% 
        layer_dense(units = nn_ls$layers[i_lyr], activation = nn_ls$actv) }
    
    # Monotonicity in output (increments; strict via softplus, else relu)
    alpha0_out <- layer_dense(object = hidden, units = 1)
    alpha_out <- layer_dense(object = hidden, units = nn_ls$p_degree, activation = nn_ls$actv_out)
    
    # Concatenate output
    output <- layer_concatenate(c(alpha0_out, alpha_out))
    
    # Define model
    model <- keras_model(inputs = c(dir_input, id_input), outputs = output)
    
    #### Estimation ####
    # Compile model
    model %>% compile(
      optimizer = custom_opt,
      loss = qt_loss
    )
    
    # Take time
    start_tm <- Sys.time()
    
    # Fit model
    history <- model %>% fit(
      x = X_train,
      y = train[i_train, bias_var],
      epochs = nn_ls$n_epochs,
      batch_size = nn_ls$n_batch,
      validation_data = list(X_valid, train[i_valid, bias_var]),
      verbose = nn_ls$nn_verbose,
      callbacks = callback_early_stopping(patience = nn_ls$n_patience,
                                          restore_best_weights = TRUE,
                                          monitor = "val_loss")
    )
    
    # Take time
    end_tm <- Sys.time()
    
    # Time needed
    runtime_est <- runtime_est + as.numeric(difftime(end_tm, start_tm, units = "mins"))
    
    # Save model
    if(model_out){ model_ls[[i_sim]] <- model }
    
    # Delete history
    rm(history)
    
    #### Prediction ####
    # Take time
    start_tm <- Sys.time()
    
    # Predict coefficients of Bernstein polynomials
    coeff_bern <- coeff_bern + predict(model, X_pred)
    
    # Take time
    end_tm <- Sys.time()
    
    # Time needed
    runtime_pred <- runtime_pred + as.numeric(difftime(end_tm, start_tm, units = "mins"))
    
    # Delete model
    rm(id_input, dir_input, station_embedding_part, hidden, 
       output, model)
    
    # Clear memory and session
    gc()
    k_clear_session()
  }
  
  # Delete data
  rm(X_train, X_valid, X_pred)
  
  # Average increments
  coeff_bern <- coeff_bern/nn_ls$n_sim
  
  # Accumulate increments
  coeff_bern <- t(apply(coeff_bern, 1, cumsum))
  
  #### Evaluation ####
  # Sum up calculated quantiles of bias (Sum of basis at quantiles times coefficients)
  q <- bern_quants(alpha = coeff_bern,
                   q_levels = q_levels)
  
  # Add up forecast and bias while ensuring positivity
  q <- pmax(X[,fc_var] - q, 0) 
  
  # Orientation has changed
  q <- t(apply(q, 1, rev))
  
  # Calculate evaluation measure of BQN forecasts
  scores_pp <- fn_scores_ens(ens = q,
                             y = X[[obs_var]],
                             skip_evals = c("e_me"),
                             scores_ens = scores_pp)
  
  # Transform ranks to uPIT
  scores_pp[["pit"]] <- scores_pp[["rank"]]/(ncol(q) + 1) - 
    runif(n = n_test,
          min = 0,
          max = 1/(ncol(q) + 1))
  scores_pp[["rank"]] <- NULL
  
  # Calculate bias of mean forecast (formula given; fc. and bias can be added up analog.)
  scores_pp[["e_me"]] <- pmax(X[,fc_var] - rowMeans(coeff_bern), 0) - X[[obs_var]]
  
  #### Output ####
  # Output list
  res <- list(f = q, 
              q_levels = q_levels,
              alpha = coeff_bern,
              nn_ls = nn_ls,
              scores_pp = scores_pp, 
              pred_vars = pred_vars,
              n_preds = length(pred_vars),
              n_train = n_train,
              n_valid = n_valid,
              n_test = n_test,
              runtime_est = runtime_est,
              runtime_pred = runtime_pred,
              center = tr_center,
              scale = tr_scale)
  
  # Include model?
  if(model_out){ res[["model_ls"]] <- model_ls }
  
  # Output
  return(res)
}