#' Breslow estimator of baseline survival
#'
#' Estimates the baseline survival function from a Cox proportional hazards model following Breslow's estimator.
#'
#' @param t (Optional) Value at which the baseline survival function is to be calculated. If \code{t = NULL} (the default), values at all unique event times is returned.
#' @param W Column name of observed predictor values (including censored opens). 
#' @param Delta Column name of censoring indicators. Note that \code{Delta = 0} is interpreted as a censored observation. 
#' @param HR Column name of hazard ratios (HR) from a Cox proportional hazards model. 
#' @param data Dataframe or named matrix containing columns \code{W}, \code{Delta}, and \code{HR}.
#'
#' @return If \code{t} is not \code{NULL}, then a scalar value for the baseline survival function at value \code{t}. Otherwise, a list with the following three elements is returned:
#' \item{W}{a vector of the unique observed failure times}
#' \item{basesurv}{baseline survival estimates}
#' \item{basehaz}{baseline cumulative hazard estimates}
#'
#' @export
#' @import magrittr
#' @importFrom dplyr group_by summarize

breslow_estimator <- function (t = NULL, W, Delta, HR, data)
{
  tj <- data[, W] # all times (censored or event)
  data %>%
    dplyr::group_by(get(W)) %>%
    dplyr::summarize(d = sum(get(Delta))) %>%
    data.frame() -> dj
  colnames(dj) <- c(W, Delta)
  dj <- dj[dj[, Delta] > 0, ] # filter to event times

  if (!is.null(t)) {
    # For times prior to first event, survival = 1
    if ((min(dj[, W]) - t) > 1E-10) {
      return(1)
    } else {
      dj <- dj[(dj[, W] - t) <= 1E-10, ] # filter to event times up to at_x
    }
  }
  tauj <- dj[, W] # event times
  atrisk_hr <- sapply(X = tauj, FUN = at_risk, t = data[, W], HR = data[, HR])
  haz0 <- dj[, Delta]/atrisk_hr
  cumhaz0 <- cumsum(haz0)
  surv0 <- exp(- cumhaz0)
  if (!is.null(t)) {
    return(surv0[length(surv0)])
  } else {
    return(list(times = tauj, basesurv = surv0, basecumhaz = cumhaz0))
  }
}