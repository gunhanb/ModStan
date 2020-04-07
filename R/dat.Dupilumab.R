#' A phase II trial of dupilumab investigating the treatment of atopic dermatitics
#'
#' A dataset containing the results from a Phase II trial with multiple dose regimens
#' examining the efficacy of dupilumab against the atopic dermatitis.
#' (ClinicalTrials.gov identifier: NCT01859988).
#'
#' @format A data frame with following coloumns
#' \describe{
#'   \item{schedule}{Schedule indicator for each patient. 1 is weekly, 2 is biweekly,
#'   and 3 is monthly dose regimen.}
#'   \item{freq}{Frequency of administration in hours for each arm}
#'   \item{dose}{Dosing amount  for each arm}
#'   \item{dose_rescaled}{Re-scaled dosing amount for each arm (for complete pooling)}
#'   \item{resp}{Percentage change from baseline in EASI score}
#'   \item{sigma}{Standard errors for each arm}
#' }
#' @source Thaci D, Simpson EL, Beck LA, et al. Efficacy and safety of dupilumab in adults with
#'  moderate-to-severe atopic dermatitis inadequately controlled by topical treatments:
#'  a randomised, placebo controlled, dose-ranging phase 2b trial. The Lancet. 2016;387:40–52.
"dat.Dupilumab"
