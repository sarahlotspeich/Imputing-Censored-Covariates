#' Breslow estimator of baseline survival
#'
#' Estimates the baseline survival function from a Cox proportional hazards model following Breslow's estimator.
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

#### SARAH: can we change that argument name data to something else? df? surv.data?
# Also, what restrictions do we have for time, delta, and hr? Is it as follows:
# time must be (stricly?) positive
# delta must be \in (0, 1)
# hr must be (strictly?) positive
# ?
breslow_estimator <- function(time, delta, hr, data) {
  # test for bad input
  if (!is.character(time)) { stop("argument time must be a character") }
  if (!is.character(delta)) { stop("argument delta must be a character") }
  if (!is.character(hr)) { stop("argument hr must be a character") }
  if (!is.data.frame(data)) { stop("argument data must be a data frame") }
  # test that data contains columns with specified names
  if (!(time %in% colnames(data))) { stop(paste("data does not have column with name ", time)) }
  if (!(delta %in% colnames(data))) { stop(paste("data does not have column with name ", delta)) }
  if (!(hr %in% colnames(data))) { stop(paste("data does not have column with name ", hr)) }
  # test for improper entries in columns of data
  #### Still deciding whether the following conditions should produce errors (stop()) or warnings
  if (any(data[, time] < 0)) { warning(paste("elements of column ", time, " must be positive")) }
  if (!all(data[, delta] %in% c(0, 1))) { warning(paste("elements of column ", delta, " must be either 0 or 1")) }
  if (any(data[, hr] < 0)) { warning(paste("elements of column ", hr, " must be inclusively between 0 and 1"))}
  
  
  #### SARAH: can you help me to add comments starting from here, or correct any of the comments I've added?
  # observed time
  tj <- data[, time]
  dj <- aggregate(formula = as.formula(paste(delta, "~", time)), FUN = sum, data = data)
  dj <- dj[dj[, delta] > 0, ]
  # observed failure times
  tauj <- dj[, time]
  riskset <- data.frame(tauj = rep(tauj, times = length(tj)),
                        tj = rep(tj, each = length(tauj)),
                        hr = rep(data[, hr], each = length(tauj)))
  riskset$atrisk <- with(riskset, tauj <= tj)
  riskset$atrisk_hr <- with(riskset, atrisk * hr)
  # denominator for baseline hazard estimate
  haz0_denom <- aggregate(formula = atrisk_hr ~ tauj, FUN = sum, data = riskset)
  
  # baseline hazard estimate
  haz0 <- dj$d / haz0_denom$atrisk_hr
  # baseline cumulative hazard estimate
  cumhaz0 <- cumsum(haz0)
  # baseline survival estimate
  surv0 <- exp(- cumhaz0)
  
  return(list(times = tauj, basesurv = surv0, basehaz = cumhaz0))
}
