#' Breslow estimator of baseline survival
#'
#' Estimates the baseline survival function from a Cox proportional hazards model following Breslow's estimator.
#'
#' @param t (Optional) Value at which the baseline survival function is to be calculated. If \code{t = NULL} (the default), values at all unique event times is returned.
#' @param W Numeric vector of observed covariate values (including censored opens). 
#' @param Delta Numeric vector of censoring indicators to accompany \code{W}. Note that \code{Delta = 0} is interpreted as a censored observation. 
#' @param HR Numeric vector of hazard ratios (HR) corresponding to values in \code{W}. 
#'
#' @return If \code{t} is not \code{NULL}, then a scalar value for the baseline survival function at value \code{t}. Otherwise, a list with the following three elements is returned:
#' \item{times}{a vector of the unique observed failure times}
#' \item{basesurv}{baseline survival estimates}
#' \item{basehaz}{baseline cumulative hazard estimates}
#'
#' @export
#' @importFrom dplyr group_by summarize

breslow_estimator <- function(t, W, Delta, HR) {
  data <- data.frame(W, Delta, HR)
  wj <- W # all times (censored or event)
  dj <- group_by(data, W)
  dj <- summarize(dj, d = sum(Delta))
  dj <- data.frame(dj)
  dj <- dj[dj$d > 0, ] # filter to event times
  
  if (!is.null(t)) {
    # For times prior to first event, survival = 1
    if ((min(dj$W) - t) > 1E-10) {
      return(1)
    } else {
      dj <- dj[(dj$W - t) <= 1E-10, ] # filter to event times up to at_x
    }
  }
  tauj <- dj$W # event times
  atrisk_hr <- sapply(X = tauj, FUN = at_risk, W = W, HR = HR)
  haz0 <- dj$d/atrisk_hr
  cumhaz0 <- cumsum(haz0)
  surv0 <- exp(- cumhaz0)
  if (!is.null(t)) {
    return(surv0[length(surv0)])
  } else {
    return(list(times = tauj, basesurv = surv0, basecumhaz = cumhaz0))
  }
}
