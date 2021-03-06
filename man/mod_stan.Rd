% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mod_stan.R
\name{mod_stan}
\alias{mod_stan}
\title{Dose-response modeling using Stan}
\usage{
mod_stan(
  dose = NULL,
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
  init = "random",
  chains = 3,
  iter = 2000,
  warmup = 1000,
  stan_seed = 1234,
  adapt_delta = 0.95
)
}
\arguments{
\item{dose}{A numerical vector specfiying the dose values}

\item{freq}{A numerical vector specfiying the frequency of administrations in hours.}

\item{schedule}{A numerical vector specfiying the schedule indicator}

\item{resp}{A numerical vector specfiying the response values}

\item{freq_ref}{A numerical value specifying the reference frequency of administration in hours}

\item{sigma}{An optional vector specfying the standard errors assosciated with `resp`.}

\item{data}{Optional data frame containing the variables given to the arguments above}

\item{N_pred}{Number of predicted dose. Default is 30.}

\item{model}{A string specifying the model used. Available options are
`CP` (complete pooling), `PP-FE` (partial pooling with fixed effects),
and `PP-RE` (partial pooling with random effects). Default is `PP-RE`.}

\item{tau_prior}{A numerical value specifying the standard deviation of the prior density
for heterogenety stdev. Default is 0.5.}

\item{tau_prior_dist}{A string specifying the prior density for the heterogeneity standard deviation,
option is `half-normal` for half-normal prior, `uniform` for uniform prior, `half-cauchy` for
half-cauchy prior.}

\item{init}{Initial values specification. The default means they are generated
randomly by `rstan`.}

\item{chains}{A positive integer specifying the number of Markov chains.
The default is 4.}

\item{iter}{A positive integer specifying the number of iterations for each chain
(including warmup). The default is 2000.}

\item{warmup}{A positive integer specifying the number of warmup (aka burnin)
iterations per chain. The default is 1000.}

\item{stan_seed}{The seed for random number of generator of Stan.}

\item{adapt_delta}{A numerical value specfying the target average proposal acceptance
probability for adaptation. See Stan manual for details. Default is 0.95. In general
you should not need to change adapt_delta unless you see a warning message about
divergent transitions, in which case you can increase adapt_delta from the
default to a value closer to 1 (e.g. from 0.95 to 0.99, or from 0.99 to 0.999, etc).}
}
\value{
an object of class `stanfit` returned by `rstan::sampling`
}
\description{
`mod_stan` fits a dose-response model with multiple schedules
 using Stan. Currently, only "emax" model is available.
}
\examples{
\dontrun{
data('dat.Dupilumab', package = "ModStan")
## Fitting an Emax model using PP-RE estimation for ED50 parameters
partial.PP-RE.Dupilumab.stan = mod_stan(dose = dose,
                                           resp = resp,
                                           schedule = schedule,
                                           freq = freq,
                                           sigma = sigma,
                                           data = dat.dupilumab,
                                           model = "PP-RE",
                                           freq_ref = 2,
                                           tau_prior_dist = "half-normal",
                                           tau_prior = 0.5,
                                           chains = 4,
                                           iter = 2000,
                                           warmup = 1000,
                                           stan.seed = 12243)
## Obtaining a small summary
print(RE.Dupilumab.stan)
## Extract the rstan fit for post-processing, eg convergence diagnostics
PP-RE.dupilumab.stanfit = partial.PP-RE.Dupilumab.stan$fit_sum

}

}
