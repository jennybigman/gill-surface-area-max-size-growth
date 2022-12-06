data {
  int<lower=1> N; //the number of observations
  int<lower=1> J; //the number of groups
  int<lower=1,upper=J> sp[N]; //vector of group indices (group identifier)
  vector[N] LogGSAcm2; //the response variable model 1
  vector[N] LogCenterMassG; // preidctor model1 (mass)
  vector[J] GP; // response variable model 2
  vector[J] CFAR; // caudal fin aspect ratio (standardized)
  matrix[J, J] d_mat; // sigma matrix
  matrix[J, J] vcov_mat; // vcov matrix


  }

 parameters {
  
  //level 1
  real global_int;
  real global_slope;
  vector[J] beta_int;
  vector[J] beta_slope;
  real<lower=0> sigma;
  
  //level 2
  real aGP; // intercept of model 2
  real bGP_slope; // slope of model 2
  real bGP_CFAR; // slope of model 2
  real<lower=0> sigma_GP; // variance model 2

  real<lower=0,upper=1>lambda_GP; // phylogenetic signal
  
  }

  transformed parameters {

  vector[N] mu_LogGSAcm2; //linear predictor model 1

  real beta_slope_ref; // intercept of first contrast
  vector[J] beta_slopes; // vector of intercepts extract from first model
  vector[J] beta_slopes_std; // vector of intercepts extract from first model standardized
  
  vector[J] mu_GP; // linear predictor model 2
  
  matrix[J, J] sigma_mat_GP;
  matrix[J, J] sigma_total_GP;
  
  
  for(n in 1:N) {

  mu_LogGSAcm2[n] = ((global_int + beta_int[sp[n]]) + (global_slope + beta_slope[sp[n]]) * LogCenterMassG[n]);
  }
  
  beta_slope_ref = global_slope;
  beta_slopes[1:J] = beta_slope_ref + beta_slope[1:J];
  beta_slopes_std = ((beta_slopes - mean(beta_slopes)) / sd(beta_slopes));

  
  mu_GP = aGP + bGP_slope * beta_slopes_std + bGP_CFAR * CFAR; // compute linear predictor for model 2
  
  sigma_mat_GP = (1-lambda_GP)*d_mat + lambda_GP*vcov_mat;
  sigma_total_GP = sigma_mat_GP * sigma_GP;

      }


model {
  
  global_int       ~ student_t(3, 0, 10);
  global_slope     ~ student_t(3, 0, 10);
  beta_int         ~ student_t(3, 0, 10);
  beta_slope       ~ student_t(3, 0, 10);
  sigma            ~ cauchy(0, 10); 
      
  LogGSAcm2 ~ normal(mu_LogGSAcm2,sigma);
  
  aGP       ~ student_t(3, 0, 10); // prior on intercept model 2
  bGP_slope   ~ student_t(3, 0, 10); // prior on slope model 2
  bGP_CFAR  ~ student_t(3, 0, 10); // prior on slope model 2
  sigma_GP  ~ cauchy(0, 10); // prior on variance model 2
  
  GP ~ normal(mu_GP, sigma_GP); 
  
  lambda_GP  ~ uniform(0,1);

  
}

generated quantities {

    real log_lik;

    log_lik = normal_lpdf(LogGSAcm2| mu_LogGSAcm2, sigma) +
              multi_normal_lpdf(GP | mu_GP, sigma_total_GP);

  }
