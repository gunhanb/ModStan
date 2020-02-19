/* 
  Burak Kursad Gunhan
  Dose-response modeling for multiple dose regimen
  Emax model
  Shared parameters (Placebo effect and Emax parameter)
  Assuming schedule-specific ED50 parameters
  Function uniform prior for ED50 parameter
*/


data {
  int<lower=1> N_obs;                           // num observations 
  int<lower=1> N_schedule;                      // num schedules
  int<lower=1> N_pred;                          // num predicted doses
  real resp[N_obs];                             // responses 
  real<lower=0> dose[N_obs];        
  real<lower=0> maxdose;
  int schedule[N_obs];                          // schedule indicator              
  real Pred_doses[N_pred];                      // predicted doses
  real<lower=0> prior_stdev_theta_0;          // prior stdev for E0
  real<lower=0> prior_stdev_theta_1;          // prior stdev for Emax
  real<lower=0> prior_stdev_sigma;          // prior stdev for sigma
}

parameters {
  real theta_0;                               // placebo effect (shared) 
  real theta_1;                               // Emax parameter (shared)
  real<lower=0, upper=1.5> theta_2_raw[N_schedule];          // ED50 parameter (not shared)
  real<lower=0> sigma;                        
}

transformed parameters{
  vector[N_obs] resp_hat;
  real<lower=0> theta_2[N_schedule];        

  
  for(i in 1:N_schedule) {
    theta_2[i] = maxdose * theta_2_raw[i];
  }
  // Dose-response: 
  for(i in 1:N_obs) 
    resp_hat[i] = theta_0 + (theta_1 * dose[i]) / (theta_2[schedule[i]] + dose[i]);
    
}

model {
  // likelihood
  for(i in 1:N_obs)  
    resp[i] ~ normal(resp_hat[i], sigma);

  // prior distributions
  sigma   ~ normal(0, prior_stdev_sigma);
  theta_0 ~ normal(0, prior_stdev_theta_0);
  theta_1 ~ normal(0, prior_stdev_theta_1);
  // Approximation to the functional uniform prior (NOT USED)
  theta_2_raw ~ lognormal(-2.5, 1.8);
}

generated quantities {
  // Predictions
  matrix[N_pred,N_schedule] resp_pred_mean; 
  vector[N_schedule] resp_pred[N_pred];
  // PP checks
  vector[N_obs] resp_rep_mean;
  vector[N_obs] resp_rep;
  // Model comparison
  vector[N_obs] log_lik;
  real theta_2_prior_raw;
  real<lower=0> theta_2_prior[N_schedule];          // ED50 parameter (not shared) <ordered vector>


  // Calculate predictions
  for(i in 1:N_pred) {
      for(j in 1:N_schedule) {
        resp_pred_mean[i,j] = theta_0 + (theta_1 * Pred_doses[i]) / (theta_2[j] + Pred_doses[i]);
        resp_pred[i,j] = normal_rng(resp_pred_mean[i,j], sigma);
      }
  }

  // Prior predictive of ED50 parameters  
  theta_2_prior_raw = lognormal_rng(-2.5, 1.8);
  theta_2_prior[2] = theta_2_prior_raw * maxdose;
  theta_2_prior[1] = theta_2_prior[2] / 2;
  theta_2_prior[3] = theta_2_prior[2] * 2;


  // PPcheck: 
  for(i in 1:N_obs) {
    resp_rep_mean[i] = theta_0 + (theta_1 * dose[i]) / (theta_2[schedule[i]] + dose[i]);
    resp_rep[i] = normal_rng(resp_rep_mean[i], sigma);
  }


  // loglikelihood contributions
  for (n in 1:N_obs)  
    log_lik[n] = normal_lpdf(resp[n]| resp_hat[n], sigma);

}



