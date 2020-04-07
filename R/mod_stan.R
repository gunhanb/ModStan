#' Dose-response modeling using Stan
#'
#' `mod_stan` fits a dose-response model with multiple schedules
#'  using Stan. Currently, only "emax" model is available.
#'
#' @export
#' @param dose A numerical vector specfiying the dose values
#' @param resp A numerical vector specfiying the response values
#' @param sigma An optional vector specfying the standard errors assosciated with `resp`.
#' @param schedule A numerical vector specfiying the schedule indicator
#' @param freq A numerical vector specfiying the frequency of administrations in hours.
#' @param N_pred Number of predicted dose. Default is 30.
#' @param model A string specifying the model used. Available options are
#' `CP` (complete pooling), `PP-FE` (partial pooling with fixed effects),
#' and `PP-RE` (partial pooling with random effects). Default is `PP-RE`.
#' @param freq_ref A numerical value specifying the reference frequency of administration in hours
#' @param data Optional data frame containing the variables given to the arguments above
#' @param tau_prior A numerical value specifying the standard deviation of the prior density
#' for heterogenety stdev. Default is 0.5.
#' @param tau_prior_dist A string specifying the prior density for the heterogeneity standard deviation,
#' option is `half-normal` for half-normal prior, `uniform` for uniform prior, `half-cauchy` for
#' half-cauchy prior.
#' @param adapt_delta A numerical value specfying the target average proposal acceptance
#' probability for adaptation. See Stan manual for details. Default is 0.95. In general
#' you should not need to change adapt_delta unless you see a warning message about
#' divergent transitions, in which case you can increase adapt_delta from the
#' default to a value closer to 1 (e.g. from 0.95 to 0.99, or from 0.99 to 0.999, etc).
#' @param iter A positive integer specifying the number of iterations for each chain
#' (including warmup). The default is 2000.
#' @param warmup A positive integer specifying the number of warmup (aka burnin)
#' iterations per chain. The default is 1000.
#' @param init Initial values specification. The default means they are generated
#' randomly by `rstan`.
#' @param chains A positive integer specifying the number of Markov chains.
#' The default is 4.
#' @param stan_seed The seed for random number of generator of Stan.
#' @return an object of class `stanfit` returned by `rstan::sampling`
#'
#' @examples
#' \dontrun{
#' data('dat.Dupilumab', package = "ModStan")
#' ## Fitting an Emax model using PP-RE estimation for ED50 parameters
#' partial.PP-RE.Dupilumab.stan = mod_stan(dose = dose,
#'                                            resp = resp,
#'                                            schedule = schedule,
#'                                            freq = freq,
#'                                            sigma = sigma,
#'                                            data = dat.dupilumab,
#'                                            model = "PP-RE",
#'                                            freq_ref = 2,
#'                                            tau_prior_dist = "half-normal",
#'                                            tau_prior = 0.5,
#'                                            chains = 4,
#'                                            iter = 2000,
#'                                            warmup = 1000,
#'                                            stan.seed = 12243)
#' ## Obtaining a small summary
#' print(RE.Dupilumab.stan)
#' ## Extract the rstan fit for post-processing, eg convergence diagnostics
#' PP-RE.dupilumab.stanfit = partial.PP-RE.Dupilumab.stan$fit_sum
#'
#' }
#'
mod_stan = function(dose = NULL,
                    freq = NULL,
                    schedule = NULL,
                    resp = NULL,
                    freq_ref = NULL,
                    sigma = NULL,
                    data = NULL,
                    N_pred = 30,
                    model = "PP-RE",
                    tau_prior = 0.5,
                    tau_prior_dist = "half-normal",
                    init = 'random',
                    chains = 3,
                    iter = 2000,
                    warmup = 1000,
                    stan_seed = 1234,
                    adapt_delta = 0.95) {
  ################ check model used
  if (model %in% c("CP", "PP-RE", "PP-FE") == FALSE) {
    stop("Function argument \"model\" must be equal to \"CP\" or \"PP-RE\" or \"PP-FE\"!!!")
  }

  ################ check reference dose regimen
  if (model == c("CP") & is.null(freq_ref) == TRUE) {
    stop("Function argument \"freq_ref\" must be provided!!!")
  }

  ################ check prior for heterogeneity parameter
  if(is.null(tau_prior_dist) == TRUE){
    stop("Function argument \"half-normal\" or \"uniform\" or \"half-cauchy\" must be specified !!!")
  }

  if(tau_prior_dist == "half-normal") { tau_prior_dist_num = 1 }
  if(tau_prior_dist == "uniform")     { tau_prior_dist_num = 2 }
  if(tau_prior_dist == "half-cauchy") { tau_prior_dist_num = 3 }

  ################ check data argument
  if (is.null(data))
    data <- sys.frame( sys.parent() )
  mf <- match.call()
  mf$data = NULL
  mf$model = NULL
  mf$tau_prior = NULL
  mf$tau_prior_dist = NULL
  mf$init = NULL
  mf$chains = NULL
  mf$iter = NULL
  mf$warmup = NULL
  mf$adapt_delta = NULL
  mf$N_pred = NULL
  mf[[1]] <- as.name("data.frame")
  mf <- eval(mf,data)
  dose     <- as.numeric(mf$dose) ## integers might overflow
  sigma     <- as.numeric(mf$sigma) ## integers might overflow
  resp <- as.numeric(mf$resp)
  freq <- as.numeric(mf$freq)
  freq <- freq[!duplicated(freq)]
  schedule <- as.numeric(mf$schedule)


  Pred_doses = seq(0, max(dose), length.out = N_pred)


  ## Create a list to be used with Stan
  ## For  PP-FE and PP-RE
  stanDat = list(N_obs = length(dose),
                 N_schedule = max(schedule),
                 resp = resp,
                 dose = dose,
                 sigma = sigma,
                 freq_ref = freq_ref,
                 schedule = schedule,
                 N_pred = N_pred,
                 Pred_doses = Pred_doses,
                 prior_stdev_theta_0 = 100,
                 prior_stdev_theta_1 = 100)




  ## Fitting the model
  if(model == "PP-RE") {
    stanDat$prior_stdev_tau = tau_prior
    stanDat$tau_prior_dist  = tau_prior_dist_num

    fit = rstan::sampling(stanmodels$PP_RE,
                          data = stanDat,
                          chains = chains,
                          seed = stan_seed,
                          iter = iter,
                          warmup = warmup,
                          control = list(adapt_delta = adapt_delta))
  }
  if(model == "PP-FE") {
    fit = rstan::sampling(stanmodels$PP_FE,
                          data = stanDat,
                          chains = chains,
                          iter = iter,
                          seed = stan_seed,
                          warmup = warmup,
                          control = list(adapt_delta = adapt_delta))
  }


  if(model == "CP") {
    Pred_doses = seq(0, max(dose), length.out = N_pred)

    stanDat_CP <- list(N_obs = length(dose),
                       resp = resp,
                       dose = dose,
                       sigma = sigma,
                       N_pred = N_pred,
                       Pred_doses = Pred_doses)

    fit = rstan::sampling(stanmodels$CP,
                          data = stanDat_CP,
                          chains = chains,
                          seed = stan_seed,
                          warmup = warmup,
                          control = list(adapt_delta = adapt_delta))
  }

  ## MODEL FINISHED
  fit_sum <- rstan::summary(fit)$summary


  Rhat.max <- max(fit_sum[,"Rhat"], na.rm = TRUE)

  if(Rhat.max > 1.1)
    warning("Maximal Rhat > 1.1. Consider increasing meta_stan MCMC parameters.")

  ## finally include a check if the Stan NuTS sample had any
  ## divergence in the sampling phase, these are not supposed to
  ## happen and can often be avoided by increasing adapt_delta
  sampler_params <- rstan::get_sampler_params(fit, inc_warmup=FALSE)
  n_divergent <- sum(sapply(sampler_params, function(x) sum(x[,'divergent__'])) )
  if(n_divergent > 0) {
    warning(paste("In total", n_divergent, "divergent transitions occured during the sampling
                  phase.\nPlease consider increasing adapt_delta closer to 1."))
  }

  out = list(fit = fit,
             fit_sum = fit_sum,
             model = model,
             data = stanDat)
  class(out) <- "mod_stan"

  return(out)

}

