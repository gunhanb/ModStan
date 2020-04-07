#' Print mod_stan object
#'
#' Takes an \code{mod_stan} object which is obtained by function \code{mod_stan} and print
#' the model and data information such as model type used in the model.
#'
#' The resulting data.frame can be used as data argument in \code{mod_stan}.
#'
#' @param x A \code{mod_stan} object.
#' @param digits An integer indicating the number of decimal places.
#' @param ... Further arguments passed to or from other methods.
#' @return The return value is invisible \code{NULL}
#' @export
print.mod_stan <- function(x, digits = 2, ...) {
  if (!is.element("mod_stan", class(x)))
    stop("Argument 'x' must be an object of class \"mod_stan\".")
  results = x$fit_sum
  cat("Dose-response modeling using ModStan\n")
  cat("E_0\n")
  print(round(results[c('E0'), -c(2, 3, 5, 7, 9, 10)], digits))
  cat("E_max\n")
  print(round(results[c('Emax'), -c(2, 3, 5, 7, 9, 10)], digits))
  cat("ED_50 (reference frequency)\n")
  print(round(results[c('ED50[2]'), -c(2, 3, 5, 7, 9, 10)], digits))

  if (x$model %in% c("PP-RE") == TRUE){
    cat("Between-schedule heterogeneity\n")
    print(round(results['tau_ED50', -c(2, 3, 5, 7, 9, 10)], digits))
  }

  return(invisible())
}
