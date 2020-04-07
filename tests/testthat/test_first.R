
context("Checking the dupilumab trial")
suppressPackageStartupMessages(library(rstan))

## Shrinkage estimation
test_that("Results are correct for fitting dose-response model using Shrinkage estimation.", {
  skip_on_cran()

  set.seed(23344)
  ## Load the dataset
  data('dat.Dupilumab', package = "ModStan")

  ## Fitting an Emax model using shrinkage estimation for ED50 parameters

  shrinkage.Dupilumab.stan  <- mod_stan(dose = dose,
                                        resp = resp,
                                        schedule = schedule,
                                        sigma = sigma,
                                        freq = freq,
                                        freq_ref = 24 * 7 * 2,
                                        data = dat.Dupilumab,
                                        model = "PP-RE",
                                        stan_seed = 2234,
                                        tau_prior =  1,
                                        tau_prior_dist = "half-normal")
  ### compare with results
  results = shrinkage.Dupilumab.stan$fit_sum

  expect_equivalent(round(results['Emax', '50%'], 2), -59.7, tolerance = 0.5)

})

## Pooling method
test_that("Results are correct for fitting dose-response model using Pooling method.", {
  skip_on_cran()

  set.seed(1111)
  ## Load the dataset
  data('dat.Dupilumab', package = "ModStan")

  ## Fitting an Emax model using Pooling
  pooling.Dupilumab.stan  <- mod_stan(dose = dose_rescaled,
                                      resp = resp,
                                      schedule = schedule,
                                      sigma = sigma,
                                      freq = freq,
                                      freq_ref = 24 * 7 * 2,
                                      data = dat.Dupilumab,
                                      model = "CP",
                                      stan_seed = 2234)
  ### compare with results
  results = pooling.Dupilumab.stan$fit_sum

  expect_equivalent(round(results['theta_0', '50%'], 2), -16.0, tolerance = 0.5)

})


## Stratified analysis
test_that("Results are correct for fitting dose-response model using Stratified analysis.", {
  skip_on_cran()

  set.seed(555)
  ## Load the dataset
  data('dat.Dupilumab', package = "ModStan")

  ## Fitting an Emax model using shrinkage estimation for ED50 parameters
  stratified.Dupilumab.stan  <- mod_stan(dose = dose,
                                         resp = resp,
                                         schedule = schedule,
                                         sigma = sigma,
                                         freq = freq,
                                         freq_ref = 24 * 7 * 2,
                                         data = dat.Dupilumab,
                                         model = "PP-FE",
                                         stan_seed = 2234)
  ### compare with results
  results = stratified.Dupilumab.stan$fit_sum

  expect_equivalent(round(results['theta_2[2]', '50%'], 2), 21.9, tolerance = 0.5)

})
