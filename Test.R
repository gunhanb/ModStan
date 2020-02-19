library(ModStan)
data("dat.Dupilumab")
head(dat.Dupilumab)

dose = dat.Dupilumab$dose
resp = dat.Dupilumab$resp
subgroup = dat.Dupilumab$subgroup
data = dat.Dupilumab
frequency = dat.Dupilumab$frequency
reference_freq = 24 * 7 * 8
N_pred = 30
model = "Pooling"
tau_prior = 0.5
tau_prior_dist = "half-normal"
chains = 4
iter = 2000
warmup = 1000
adapt_delta = 0.95
Emax_pooling.model           <- stan_model(file = "inst/stan/Pooling.stan")
Emax_shrinkage.model           <- stan_model(file = "inst/stan/Shrinkage.stan")

data('dat.Dupilumab', package = "ModStan")
 ## Fitting an Emax model using shrinkage estimation for ED50 parameters
 shrinkage.Dupilumab.stan  <- mod_stan(dose = dose,
                                        resp = resp,
                                        subgroup = subgroup,
                                        frequency = frequency,
                                        reference_freq = 24 * 7 * 2,
                                        data = dat.Dupilumab,
                                        model = "Stratified",
                                        tau_prior =  1)
 ## Obatining a small summary
 print(shrinkage.Dupilumab.stan)
 ## Extract the rstan fit for post-processing, eg convergence diagnostics
 shrinkage.dupilumab.stanfit = shrinkage.Dupilumab.stan$fit



mod_stan = function(dose = NULL,
                    frequency = NULL,
                    subgroup = NULL,
                    resp = NULL,
                    reference_freq = NULL,
                    data = NULL,
                    N_pred = 30,
                    model = "Shrinkage",
                    tau_prior = 0.5,
                    tau_prior_dist = "half-normal",
                    chains = 4,
                    iter = 2000,
                    warmup = 1000,
                    adapt_delta = 0.95) {
  ################ check model used
  if (model %in% c("Pooling", "Shrinkage", "Stratified") == FALSE) {
    stop("Function argument \"model\" must be equal to \"Pooling\" or \"Shrinkage\" or \"Stratified\"!!!")
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


  ### Initial values
  inits = list(list("theta_0" = -20, "theta_1" = -60),
               list("theta_0" = -25, "theta_1" = -70),
               list("theta_0" = -35, "theta_1" = -80))


  ## Fitting the model
  if(model == "Shrinkage") {
    ##For Shrinkage
    stanDat$prior_stdev_tau = tau_prior
    stanDat$tau_prior_dist = tau_prior_dist_num

    fit = rstan::sampling(Emax_shrinkage.model,
                          data = stanDat,
                          chains = chains,
                          iter = iter,
                          init = inits,
                          warmup = warmup,
                          control = list(adapt_delta = adapt_delta))
  }
  if(model == "Stratified") {
    fit = rstan::sampling(stanmodels$Stratified,
                          data = stanDat,
                          chains = chains,
                          iter = iter,
                          warmup = warmup,
                          control = list(adapt_delta = adapt_delta))
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

    fit = rstan::sampling(Emax_pooling.model,
                          data = stanDat_pooling,
                          chains = chains,
                          iter = iter,
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






