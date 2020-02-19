/*
  Burak Kursad Gunhan
  Dose-response modeling
  Emax model
  Using functional uniform priors
*/
data {
  int<lower=1> N_obs;                           // num observations
  int<lower=1> N_pred;                           // num predicted doses
  real resp[N_obs];                             // responses
  real<lower=0> dose[N_obs];
  real<lower=0> maxdose;
  real Pred_doses[N_pred];                         // predicted doses
  real<lower=0> prior_stdev_sigma;                 // prior stdev for scale parameter
}

parameters {
  real theta_0;
  real theta_1;                               // Emax parameter
  real<lower=0,upper=1.5> theta_2_raw;                      // ED50 parameter
  real<lower=0> sigma;
}

transformed parameters{
  vector[N_obs] resp_hat;
  real<lower=0> theta_2;          // ED50 parameter (not shared) <ordered vector>


  theta_2 = (maxdose * theta_2_raw);

  // Dose-response: Emax model
  for(i in 1:N_obs)
    resp_hat[i] = theta_0 + (theta_1 * dose[i])/(theta_2 + dose[i]);


}

model {
  // prior distributions
  sigma   ~ normal(0, prior_stdev_sigma);
  theta_0 ~ normal(0, 100);
  theta_1 ~ normal(0, 100);
  theta_2_raw ~ lognormal(-2.5, 1.8);

  // likelihood
  resp ~ normal(resp_hat, sigma);
  // Incorporating the functional uniform prior for ED50
}


generated quantities {
  vector[N_pred] resp_pred_mean; // Predicted responses
  real resp_pred[N_pred];                             // responses
  vector[N_obs] log_lik;        // loglikelihood contribution
  real<lower=0> theta_2_new[3];          // ED50 parameter (not shared) <ordered vector>

  // Different dose regimens
  theta_2_new[1] = theta_2 / 2;
  theta_2_new[2] = theta_2;
  theta_2_new[3] = theta_2 * 2;


  // Calculate predictions
  for(i in 1:N_pred) {
    resp_pred_mean[i] = theta_0 + (theta_1 * Pred_doses[i])/(theta_2 + Pred_doses[i]);
    resp_pred[i] = normal_rng(resp_pred_mean[i], sigma);
  }


  // For model comparison
  for (n in 1:N_obs)
    log_lik[n] = normal_lpdf(resp[n]| resp_hat[n], sigma);

}

