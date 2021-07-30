#' Breslow estimator of baseline survival
#'
#' Estimates the baseling survival function from a Cox proportional hazards model following Breslow's esimator.
#'
#' @param time String column name for the observed times.
#' @param delta String column name for the censoring indicator of the covariate.
#' @param hr String column name for the hazard ratios calculated from a \code{coxph} model fit.
#' @param data Datafrane containing columns \code{time}, \code{delta}, and \code{hr}.
#'
#' @return A dataframe of the same length and in the same order as \code{data} with the following columns:
#' \item{times}{unique values of \code{data[, time]}}
#' \item{basesurv}{baseline survival estimates}
#' \item{basehaz}{baseline cumulative hazard estimates}
#'
#' @export

breslow_estimator <- function(time, delta, hr, data) {
  tj <- data[, time]
  dj <- aggregate(formula = as.formula(paste(delta, "~", time)), FUN = sum, data = data)
  dj <- dj[dj[, delta] > 0, ]
  tauj <- dj[, time]
  riskset <- data.frame(tauj = rep(tauj, times = length(tj)),
                        tj = rep(tj, each = length(tauj)),
                        hr = rep(data[, hr], each = length(tauj)))
  riskset$atrisk <- with(riskset, tauj <= tj)
  riskset$atrisk_hr <- with(riskset, atrisk * hr)
  haz0_denom <- aggregate(formula = atrisk_hr ~ tauj, FUN = sum, data = riskset)
  haz0 <- dj$d / haz0_denom$atrisk_hr
  cumhaz0 <- cumsum(haz0)
  surv0 <- exp(- cumhaz0)
  return(list(times = tauj, basesurv = surv0, basehaz = cumhaz0))
}
