library(rstan)
library(dplyr)


simulate_BASH_LMM_data <- function(min_obs_per_indiv,
                                   max_obs_per_indiv,
                                   prob_confounder1,
                                   cra_rate,
                                   msl_rate,
                                   N,
                                   beta_0,
                                   beta_cra,
                                   beta_msl,
                                   beta_time,
                                   beta_US,
                                   sd_noise,
                                   cra_delay,
                                   msl_delay,
                                   cra_retain_rate, 
                                   msl_retain_rate,
                                   theta_scale_prior,
                                   max_lag,
                                   individual_sd,
                                   output_dir){
  
  # the total number of observations per individual
  n_i = sample(min_obs_per_indiv:max_obs_per_indiv,N,replace = T)
  
  indiv = c()
  id = 1
  for(ni in n_i){
    indiv = append(indiv, rep(id,ni))
    id = id+1
  }
  
  #site ID and days_since_activation vector
  id = c()
  days_since_activation = c()
  confounder1 = c()
  for(i in 1:N){
    id = append(id, rep(i,n_i[i]))
    days_since_activation = append(days_since_activation, 1:n_i[i])
    if (rnorm(1) < qnorm(prob_confounder1)){
      confounder1 = append(confounder1, rep(1,n_i[i]))
    }else{
      confounder1 = append(confounder1, rep(0,n_i[i]))
    }
  }
  
  print("Simulating engagements")
  
  
  #engagement flag variables
  x_cra = sample(0:1,sum(n_i),replace=TRUE,prob=c(1-cra_rate,cra_rate))
  x_msl = sample(0:1,sum(n_i),replace=TRUE,prob=c(1-msl_rate,msl_rate))
  #input data frame
  X = data.frame(SITE_NUM = id, time=days_since_activation, ecdp_eng_clinical_flag = x_cra, ecdp_eng_medical_flag = x_msl,US=confounder1)
  
  #simulate output 
  
  #random intercept per individual
  z_0 = rnorm(N,mean = 0, sd = individual_sd)
  
  #STAN Input data
  N = nrow(X)
  N_indivs = length(unique(X$SITE_NUM))
  indiv = X$SITE_NUM
  
  
  control_vars = c("time", "US")
  
  # Create a control variables data frame
  X_ctrl = X %>%
    select(time, US) %>% 
    mutate(time = log(time+1))
  
  
  drop_names = c("SITE_NUM")
  drop_names =append(drop_names,names(X_ctrl))
  
  # Identify media variables
  media_vars <- setdiff(names(X),drop_names)
  num_media = length(media_vars)
  
  #2d media array with grouping var
  X_eng = X[,append("SITE_NUM",media_vars)]
  
  # Setup 3D array for media variables with lags
  X_media <- array(0, c(nrow(X_eng), num_media, max_lag))
  
  # Fill the unlagged time series for media
  X_media[, , 1] <- as.matrix(X_eng[, media_vars])
  
  # Compute lagged values and store in the 3D array
  for (lag in 1:(max_lag-1)){
    for (media in media_vars){
      X_eng = X_eng %>% group_by(SITE_NUM) %>% 
        mutate("{media}_lag_{lag}" := dplyr::lag(!!sym(media),n=lag)) %>% 
        replace(is.na(.), 0) 
    }
    X_media[,,lag+1] = as.array(as.matrix(X_eng[,(ncol(X_eng)-num_media+1):(ncol(X_eng))]))
  }
  
  # Initialize vectors
  y <- numeric(N)
  epsilon <- rnorm(N,sd = sd_noise)
  delay <- c(cra_delay, msl_delay)
  retain_rate <- c(cra_retain_rate, msl_retain_rate)
  
  # Initialize matrices for cumulative effects with Hill transformation
  cum_effects_hill <- matrix(0, nrow = N, ncol = num_media)
  
  # Placeholder coefficients for media and control variables (define appropriately)
  beta_medias <- c(beta_cra, beta_msl)
  gamma_ctrl <- c(beta_time, beta_US)
  
  
  # Hill function
  Hill <- function(t, ec, slope) {
    return(1 / (1 + (t / ec)^(-slope)))
  }
  
  # Adstock transformation with a vector of weights
  Adstock <- function(t, weights) {
    return(sum(t * weights))
  }
  
  
  print("Simulating response variable")
  # simulate response variable (normal assumption clipped > 0) with specified parameters
  for (nn in 1:N) {
    for (media in 1:num_media) {
      lag_weights <- numeric(max_lag)
      for (lag in 1:max_lag) {
        lag_weights[lag] <- retain_rate[media]^((lag - 1 - delay[media])^2)
      }
      cum_effect <- Adstock(X_media[nn, media, ], lag_weights)
      cum_effects_hill[nn, media] <- cum_effect
    }
    y[nn] <- rpois(1,exp(beta_0 +
                           sum(cum_effects_hill[nn, ] * beta_medias) +
                           sum(X_ctrl[nn,] * gamma_ctrl) + z_0[indiv[nn]] + epsilon[nn]  ) )
  }
  
  write.csv(y, paste0(output_dir,"/y.csv"))
  write.csv(X_ctrl, paste0(output_dir,"/X_ctrl.csv"))
  write.csv(X_eng, paste0(output_dir,"/X_eng.csv"))
  write.csv(cum_effects_hill, paste0(output_dir,"/cum_effects_hill.csv"))
  write.csv(z_0[indiv], paste0(output_dir,"/indiv_effects.csv"))
  write.csv(epsilon, paste0(output_dir,"/epsilon.csv"))
  
  return(list("Y"=y,
              "cum_effects_hill"=cum_effects_hill,
              "N"=N,
              "max_lag"=max_lag,
              "num_media"=num_media,
              "X_media"=X_media,
              "X_eng"=X_eng,
              "num_ctrl"=ncol(X_ctrl),
              "X_ctrl"=X_ctrl,
              "beta_medias"=beta_medias,
              "gamma_ctrl"=gamma_ctrl,
              "retain_rate"=retain_rate,
              "delay"=delay,
              "indiv"=indiv,
              "N_indivs"=N_indivs
  ))
}


input_params <- list(min_obs_per_indiv=15,
                     max_obs_per_indiv=400,
                     prob_confounder1=0.5,
                     cra_rate=0.2,
                     msl_rate=0.1,
                     N=100,
                     beta_0=-1.5,
                     beta_cra=0.4,
                     beta_msl=0.2,
                     beta_time=-0.14,
                     beta_US=-0.33,
                     sd_noise=0.4,
                     cra_delay=1,
                     msl_delay=3,
                     cra_retain_rate=0.5, 
                     msl_retain_rate=0.8,
                     priormean_icept=-1.5,
                     theta_scale_prior=0.5,
                     max_lag=13,
                     individual_sd=0.7,
                     n_iter=2000,
                     seed=9,
                     input_stan="/BASH_Poisson_LMM_no_hill.stan",
                     output_dir = "simulation_results")
