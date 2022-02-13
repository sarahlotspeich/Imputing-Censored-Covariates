#' Calculate size of risk set at \code{t} based on a vector of times \code{W}
#'
#' Calculates the size of the risk set at specified times. Alternatively, if \code{HR} is not \code{NULL} then it sums up HR of at-risk people (for use with \code{breslow_estimator()}).
#'
#' @param t Value at which the size of the risk set is to be calculated. 
#' @param W Numeric vector of observed predictor values (including censored opens). 
#' @param HR (Optional) Numeric vector of hazard ratios (HR) corresponding to values in \code{W}. 
#'
#' @return A scalar
#'
#' @export
#' 
at_risk <- function(t, W, HR = NULL) {
  if (!is.null(HR)) {
    sum((W >= t) * HR)
  } else {
    sum(W >= t)
  }
}