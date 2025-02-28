//poisson_log_adstock_mixed_model.stan
functions {
  
  // the Hill function
  real Hill(real t, real ec, real slope) {
    return 1 / (1 + (t / ec)^(-slope));
  }
  
  // the adstock transformation with a vector of weights
  real unweighted_adstock(row_vector t, row_vector weights) {
    
    return dot_product(t, weights);
  }
}





// mixed model with Poisson likelihood
data {
  
  // number of observations
  int < lower = 1 > N; 
  
  // number of non-engagement columns in X matrix
  int < lower = 1 > num_ctrl; 
  
  // number of individuals
  int < lower = 1 > N_indivs; 
  
  // subscripts indexing individuals
  int < lower = 1 > indiv[N]; 
  
  // count of events
  int < lower=0 > y[N]; 
  
  // X matrix without ones for intercept in column 1
  matrix[N, num_ctrl] X_ctrl; 
  
  // the maximum duration of lag effect, in weeks
  int<lower=1> max_lag; 
  // the number of media channels
  int<lower=1> num_media; 
  
  // 3D array of media variables
  //    X: time [1:N]
  //    Y: engagement type [1:num_media]
  //    Z: lagged time [1:max_lag]
  
  row_vector[max_lag] X_media[N, num_media];
  
}



transformed data {
  
  // QR reparameterization for X_ctrl fixed matrix
  matrix[N, num_ctrl] Q_ast; 
  matrix[num_ctrl, num_ctrl] R_ast;
  matrix[num_ctrl, num_ctrl] R_ast_inverse;
  
  // thin and scale the QR decomposition
  Q_ast = qr_Q(X_ctrl)[, 1:num_ctrl] * sqrt(N - 1);
  R_ast = qr_R(X_ctrl)[1:num_ctrl, ] / sqrt(N - 1);
  R_ast_inverse = inverse(R_ast);
  
}



parameters {
  
  real alpha; // coefficient of variation of random effects
  real < lower=0 > sigma_indiv; // SD of random effects
  vector[num_ctrl] theta;      // coefficients on Q_ast
  vector[N_indivs] xi; // random effects with standard normal prior
  
  // the coefficients for media variables
  vector<lower=0>[num_media] beta_medias;
  
  // the retention rate and delay parameter for the adstock transformation of
  // each media
  vector<lower=0,upper=1>[num_media] retain_rate;
  vector<lower=0,upper=max_lag-1>[num_media] delay;
  
  // // Hill parameters
  // vector<lower=0,upper=1>[num_media] ec;
  // vector<lower=0>[num_media] slope;
}

transformed parameters {
  
  vector[N] eta = Q_ast * theta;
  vector[N_indivs] w = alpha + xi * sigma_indiv;
  matrix[num_media, max_lag] lag_weights; // precomputed lag weights
  row_vector[num_media] cum_effects[N];
  row_vector[num_media] cum_effects_hill[N];
  
  // Precompute lag weights
  for (media in 1:num_media) {
    for (lag in 1:max_lag) {
      
      lag_weights[media, lag] = pow(retain_rate[media], 
                                    (lag - 1 - delay[media]) ^ 2);
    }
  }
  // Compute cumulative media effects
  for (nn in 1:N) {
    for (media in 1:num_media) {
      cum_effects_hill[nn, media] = unweighted_adstock(X_media[nn,media],             
                                                       lag_weights[media]);
      // cum_effects_hill[nn, media] = Hill(cum_effects[nn, media], ec[media], slope[media]);
    }
  }
}



model {
  
  // half-cauchy prior, prior values taken from DIAL models
  sigma_indiv ~ cauchy(0.7, 0.4); 
  alpha ~ normal(-3, 1);
  theta ~ normal(0,1);
  xi ~ normal(0, 1);
  // slope ~ normal(1, 1);
  // ec ~ beta(2,2);
  beta_medias ~ normal(0, 1);
  retain_rate ~ beta(4,2);
  delay ~ uniform(0, max_lag);
  
  
  
  for(i in 1:N) {
    
    y[i] ~ poisson_log(eta[i] + 
                         dot_product(cum_effects_hill[i,], beta_medias)+ w[indiv[i]]);
  }
}



generated quantities {
  
  vector[num_ctrl] beta;
  beta = R_ast_inverse * theta; // coefficients on x
  
}