/*
  Burak Kursad Gunhan
  Dose-response modeling for multiple dose regimen
  Shared parameters (Placebo effect and Eaqmx parameter)
  Assuming schedule-specific random effects ED50 parameters
*/


data {
  int<lower=1> N_obs;                           // num of observations
  int<lower=1> N_schedule;                      // num of schedules
  int<lower=1> N_pred;                          // num of predicted doses
  real resp[N_obs];                             // responses
  real<lower=0> dose[N_obs];                    // doses
  real<lower=0> sigma[N_obs];                    // standard errors
  int schedule[N_obs];                          // schedule indicator
  int<lower=0> freq_ref;
  real<lower=0> freq[N_schedule];                    // doses
  real Pred_doses[N_pred];                      // predicted doses
}
parameters {
  real E0;                                      // placebo effect (shared)
  real Emax;                                    // Emax parameter (shared)
  real log_ED50_raw[N_schedule];                // re-scaled log(ED50) parameters
  real<lower=0, upper=1.5> mu_ED50_raw;                    // mean of log(ED50) random effects
  real<lower=0> tau_ED50;                       // between-schedule heteroegeneity
}
transformed parameters{
  real mu_ED50;
  real log_ED50[N_schedule];
  real<lower=0> ED50[N_schedule];
  vector[N_obs] resp_hat;

  mu_ED50 = log(mu_ED50_raw * max(dose));
  for(i in 1:N_schedule)
    log_ED50[i] = mu_ED50 + log_ED50_raw[i] * tau_ED50;
  // Taking exponentials and rescaling ED50 parameters
  for(i in 1:N_schedule)
    ED50[i] = exp(log_ED50[i]) * (freq[i]/ freq_ref);

  // Dose-response: Emax model
  for(i in 1:N_obs)
    resp_hat[i] = E0 + (Emax * dose[i]) / (ED50[schedule[i]] + dose[i]);
}
model {
  // random effects
  log_ED50_raw ~ normal(0, 1);  // implies log(ED50) ~ normal(mu_ED50, tau_ED50)
  // likelihood
  resp ~ normal(resp_hat, sigma);
  // prior distributions
  E0 ~ normal(0, 100);
  Emax ~ normal(0, 100);
  // approximation to the functional uniform prior
  mu_ED50_raw  ~ lognormal(-2.5, 1.8);
  tau_ED50  ~ normal(0, 1);
}

generated quantities {
  matrix[N_pred,N_schedule] resp_pred_mean;           // Predictions
  vector[N_obs] log_lik;                              // Model comparison

  // Calculate predictions
  for(i in 1:N_pred) {
    for(j in 1:N_schedule) {
      resp_pred_mean[i,j] = E0 + (Emax * Pred_doses[i])/(ED50[j] + Pred_doses[i]);
    }
  }

 // loglikelihood contributions
  for (n in 1:N_obs)
    log_lik[n] = normal_lpdf(resp[n]| resp_hat[n], sigma[n]);

}
