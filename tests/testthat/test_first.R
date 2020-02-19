
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
                                        subgroup = subgroup,
                                        frequency = frequency,
                                        data = dat.Dupilumab,
                                        model = "Shrinkage",
                                        beta_prior =  1,
                                        beta_prior_dist = "half-normal")
  ### compare with results
  results = shrinkage.Dupilumab.stan$fit_sum

  expect_equivalent(round(results['theta_1', '50%'], 2), -69.3, tolerance = 0.5)

})

## Pooling method
test_that("Results are correct for fitting dose-response model using Pooling method.", {
  skip_on_cran()

  set.seed(1111)
  ## Load the dataset
  data('dat.Dupilumab', package = "ModStan")

  ## Fitting an Emax model using Pooling
  pooling.Dupilumab.stan  <- mod_stan(dose = dose,
                                        resp = resp,
                                        subgroup = subgroup,
                                        frequency = frequency,
                                        data = dat.Dupilumab,
                                        reference_freq = 24 * 7 * 2,
                                        model = "Pooling",
                                        beta_prior =  1,
                                        beta_prior_dist = "half-normal")
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
                                        subgroup = subgroup,
                                        frequency = frequency,
                                        data = dat.Dupilumab,
                                        model = "Stratified",
                                        beta_prior =  1,
                                        beta_prior_dist = "half-normal")
  ### compare with results
  results = stratified.Dupilumab.stan$fit_sum

  expect_equivalent(round(results['theta_2[2]', '50%'], 2), 21.9, tolerance = 0.5)

})
