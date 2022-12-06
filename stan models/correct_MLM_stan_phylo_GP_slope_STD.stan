  data {
      int<lower=1> N; //the number of observations
      int<lower=1> J; //the number of groups
      int<lower=1,upper=J> sp[N]; //vector of group indices (group identifier)
      vector[N] LogGSAcm2; //the response variable model 1
      vector[N] LogCenterMassG; // preidctor model1 (mass)
      vector[J] GP_slope; // response variable model 2
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
       real aGP_slope; // intercept of model 2
       real bGP_slope; // slope of model 2
       real<lower=0> sigma_GP_slope; // variance model 2
      
       real<lower=0,upper=1>lambda_GP_slope; // phylogenetic signal

      }

  transformed parameters {

      vector[N] mu_LogGSAcm2; //linear predictor model 1

      real beta_slope_ref; // intercept of first contrast
      vector[J] beta_slopes; // vector of intercepts extract from first model
      vector[J] beta_slopes_std; // vector of intercepts extract from first model standardized
      
      vector[J] muGP_slope; // linear predictor model 2
  
      matrix[J, J] sigma_mat_GP_slope;
      matrix[J, J] sigma_total_GP_slope;
      
      for(n in 1:N) {
      
      mu_LogGSAcm2[n] = ((global_int + beta_int[sp[n]]) + (global_slope + beta_slope[sp[n]]) * LogCenterMassG[n]);
      }  
        
      beta_slope_ref = global_slope;
      beta_slopes[1:32] = beta_slope_ref + beta_slope[1:32];
      beta_slopes_std = ((beta_slopes - mean(beta_slopes)) / sd(beta_slopes));
    
      muGP_slope = aGP_slope + bGP_slope * beta_slopes_std; // compute linear predictor for model 2

      sigma_mat_GP_slope = (1-lambda_GP_slope)*d_mat + lambda_GP_slope*vcov_mat;
      sigma_total_GP_slope = sigma_mat_GP_slope * sigma_GP_slope;
      
      }

  model {

      global_int       ~ student_t(3, 0, 10);
      global_slope     ~ student_t(3, 0, 10);
      beta_int         ~ student_t(3, 0, 10);
      beta_slope       ~ student_t(3, 0, 10);
      sigma            ~ cauchy(0, 10); 
          
      LogGSAcm2 ~ normal(mu_LogGSAcm2,sigma);
      
      aGP_slope      ~ student_t(3, 0, 10); // prior on intercept model 2
      bGP_slope      ~ student_t(3, 0, 10); // prior on slope model 2
      sigma_GP_slope  ~ cauchy(0, 10); // prior on variance model 2
      
      GP_slope ~ normal(muGP_slope, sigma_GP_slope); 

      lambda_GP_slope  ~ uniform(0,1);

      }


  generated quantities {

      real log_lik;

      log_lik =
      normal_lpdf(LogGSAcm2 | mu_LogGSAcm2, sigma) +
      multi_normal_lpdf(GP_slope | muGP_slope, sigma_total_GP_slope);
      
      }
