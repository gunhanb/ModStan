
context("Checking the dupilumab trial")
suppressPackageStartupMessages(library(rstan))

## Fitting a binomial-normal hierachical model with Smith parametrization
test_that("Results are correct for fitting dose-response model suing Shrinkage estimation.", {
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
                                        beta_prior =  1)
  ### compare with results
  results = shrinkage.Dupilumab.stan$fit_sum

  expect_equivalent(round(results['theta_1', '50%'], 2), -69.3, tolerance = 0.5)

})


