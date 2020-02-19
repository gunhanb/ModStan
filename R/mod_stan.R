#' Dose-response modeling using Stan
#'
#' `mod_stan` fits a dose-response model with subgroups using Stan. Currently, only "emax"
#' model is available.
#'
#' @export
#' @param dose A numerical vector specfiying the dose values
#' @param resp A numerical vector specfiying the response values
#' @param subgroup A numerical vector specfiying the subgroup indicator
#' @param frequency A numerical vector specfiying the frequency of adminsitration in hours.
#' This is only needed when subgroups refer to dose regimens.
#' @param N_pred Number of predicted dose. Default is 30.
#' @param reference_freq A numerical value specifying the frequency of administration in hours
#' for the reference dose regimen. Needed for `Pooling` method.
#' @param data Optional data frame containing the variables given to the arguments above
#' @param beta_prior A numerical value specifying the standard deviation of the prior density
#' for heterogenety stdev. Default is 0.5.
#' @param beta_prior_dist A string specifying the prior density for the heterogeneity standard deviation,
#' option is `half-normal` for half-normal prior, `uniform` for uniform prior, `half-cauchy` for
#' half-cauchy prior.
#' @param model A string specifying the model used. Available options are `Pooling`,
#'  `Stratified`, and `Shrinkage`. Default is `Shrinkage`.
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
#' @return an object of class `stanfit` returned by `rstan::sampling`
#'
#' @examples
#' \donttest{
#' data('dat.Dupilumab', package = "ModStan")
#' ## Fitting an Emax model using shrinkage estimation for ED50 parameters
#' partial.shrinkage.Dupilumab.stan = mod_stan(dose = dose,
#'                                            resp = resp,
#'                                            subgroup = subgroup,
#'                                            frequency = frequency,
#'                                            data = dat.dupilumab,
#'                                            model = "Shrinkage",
#'                                            reference_freq = NULL,
#'                                            beta_prior_dist = "half-normal",
#'                                            beta_prior = 0.5,
#'                                            chains = 4,
#'                                            iter = 2000,
#'                                            warmup = 1000)
#' ## Obtaining a small summary
#' print(shrinkage.Dupilumab.stan$fit, pars = c("theta_0", "theta_1", "theta_2", "beta_theta_2", "sigma"))
#' ## Extract the rstan fit for post-processing, eg convergence diagnostics
#' shrinkage.dupilumab.stanfit = partial.shrinkage.Dupilumab.stan$fit_sum
#'
#' }
#'
mod_stan = function(dose = NULL,
                    frequency = NULL,
                    subgroup = NULL,
                    resp = NULL,
                    reference_freq = NULL,
                    data = NULL,
                    N_pred = 30,
                    model = "Shrinkage",
                    beta_prior = 0.5,
                    beta_prior_dist = "half-normal",
                    init = 'random',
                    chains = 3,
                    iter = 2000,
                    warmup = 1000,
                    adapt_delta = 0.95) {
  ################ check model used
  if (model %in% c("Pooling", "Shrinkage", "Stratified") == FALSE) {
    stop("Function argument \"model\" must be equal to \"Pooling\" or \"Shrinkage\" or \"Stratified\"!!!")
  }

  ################ check reference dose regimen
  if (model == c("Pooling") & is.null(reference_freq) == TRUE) {
    stop("Function argument \"reference_freq\" must be provided!!!")
  }

  ################ check prior for heterogeneity parameter
  if(is.null(beta_prior_dist) == TRUE){
    stop("Function argument \"half-normal\" or \"uniform\" or \"half-cauchy\" must be specified !!!")
  }

  if(beta_prior_dist == "half-normal") { beta_prior_dist_num = 1 }
  if(beta_prior_dist == "uniform")     { beta_prior_dist_num = 2 }
  if(beta_prior_dist == "half-cauchy") { beta_prior_dist_num = 3 }

  ################ check data argument
  if (is.null(data))
    data <- sys.frame( sys.parent() )
  mf <- match.call()
  mf$data = NULL
  mf$model = NULL
  mf$beta_prior = NULL
  mf$beta_prior_dist = NULL
  mf$init = NULL
  mf$chains = NULL
  mf$iter = NULL
  mf$warmup = NULL
  mf$adapt_delta = NULL
  mf$N_pred = NULL
  mf[[1]] <- as.name("data.frame")
  mf <- eval(mf,data)
  dose     <- as.numeric(mf$dose) ## integers might overflow
  resp <- as.numeric(mf$resp)
  frequency <- as.numeric(mf$frequency)
  subgroup <- as.numeric(mf$subgroup)


  Pred_doses = seq(0, max(dose), length.out = N_pred)


  ## Create a list to be used with Stan
  ## For  stratified
  stanDat = list(N_obs = length(dose),
                 N_schedule = max(subgroup),
                 resp = resp,
                 dose = dose,
                 schedule = subgroup,
                 N_pred = N_pred,
                 maxdose = max(dose),
                 Pred_doses = Pred_doses,
                 prior_stdev_theta_0 = 100,
                 prior_stdev_theta_1 = 100,
                 prior_stdev_sigma = 100)




  ## Fitting the model
  if(model == "Shrinkage") {
    ##For Shrinkage
    stanDat$prior_stdev_beta = beta_prior
    stanDat$beta_prior_dist = beta_prior_dist_num

    fit = rstan::sampling(stanmodels$Shrinkage,
                          data = stanDat,
                          chains = chains,
                          iter = iter,
                          init = init,
                          warmup = warmup,
                          control = list(adapt_delta = adapt_delta))
  }
  if(model == "Stratified") {
    fit = rstan::sampling(stanmodels$Stratified,
                          data = stanDat,
                          chains = chains,
                          iter = iter,
                          warmup = warmup,
                          control = list(adapt_delta = adapt_delta),
                          init = init)
  }


  if(model == "Pooling") {
    ## For pooling
    for(i in 1:max(subgroup)) {
      data$dose[subgroup == i] = data$dose[subgroup == i] * (reference_freq / data$frequency[subgroup == i][1])
    }

    Pred_doses = seq(0, max(dose), length.out = N_pred)

    stanDat_pooling <- list(N_obs = length(dose),
                            resp = resp,
                            dose = dose,
                            N_pred = N_pred,
                            maxdose = max(dose),
                            Pred_doses = Pred_doses,
                            prior_stdev_sigma = 100)

    fit = rstan::sampling(stanmodels$Pooling,
                          data = stanDat_pooling,
                          chains = chains,
                          iter = iter,
                          warmup = warmup,
                          control = list(adapt_delta = adapt_delta),
                          init = init)
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

