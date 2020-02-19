/*
  Burak Kursad Gunhan
  Dose-response modeling for multiple dose regimen
  Emax model
  Shared parameters (Placebo effect)
  Assuming schedule-specific ED50 parameters
  Putting random effects on Emax parameter
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
  real<lower=0> prior_stdev_beta;          // prior stdev for sigma
}

parameters {
  real theta_0;                               // placebo effect (shared)
  real theta_1;
  real theta_2_raw;                               // Emax parameter (shared)
  real<lower=0> sigma;
  real<lower=0> mu_theta_2_raw;
  real<lower=0> beta_theta_2;
}

transformed parameters{
  real log_mu_theta_2;
  vector[N_obs] resp_hat;
  real log_theta_2[N_schedule];          // ED50 parameter (not shared) <ordered vector>
  real<lower=0> theta_2[N_schedule];          // ED50 parameter (not shared) <ordered vector>

  log_mu_theta_2 = log(mu_theta_2_raw * maxdose);

  // Reference dose regimen
  log_theta_2[2] = log_mu_theta_2;
  // Borrowing strength from other dose regimens
  log_theta_2[1] = (log_mu_theta_2 + theta_2_raw * beta_theta_2);
  log_theta_2[3] = (log_mu_theta_2 + theta_2_raw * beta_theta_2);

  // Back transformation
  theta_2[2] = exp(log_theta_2[1]);
  theta_2[1] = exp(log_theta_2[2]) / 2;
  theta_2[3] = exp(log_theta_2[2]) * 2;

  // Dose-response: Emax model
  for(i in 1:N_obs)
    resp_hat[i] = theta_0 + (theta_1 * dose[i]) / (theta_2[schedule[i]] + dose[i]);

}

model {
  // random effects
  theta_2_raw ~ normal(0, 1);  // implies theta_2 ~ normal(mu_theta_2, beta_theta_2)

  // likelihood
  resp ~ normal(resp_hat, sigma);

  // prior distributions
  sigma   ~ normal(0, prior_stdev_sigma);
  theta_0 ~ normal(0, prior_stdev_theta_0);
  theta_1 ~ normal(0, prior_stdev_theta_1);
  // approximation to the function uniform prior
  mu_theta_2_raw  ~ lognormal(-2.5, 1.8);
  beta_theta_2  ~ normal(0, prior_stdev_beta);

}
generated quantities {
  // Predictions
  matrix[N_pred,N_schedule] resp_pred_mean;
  // Model comparison
  vector[N_obs] log_lik;


  // Calculate predictions
  for(i in 1:N_pred) {
    for(j in 1:N_schedule) {
      resp_pred_mean[i,j] = theta_0 + (theta_1 * Pred_doses[i])/(theta_2[j] + Pred_doses[i]);
      }
    }

 // loglikelihood contributions
  for (n in 1:N_obs)
    log_lik[n] = normal_lpdf(resp[n]| resp_hat[n], sigma);

}

