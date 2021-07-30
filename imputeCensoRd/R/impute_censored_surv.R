#' Impute survival at censored time
#'
#' Imputes survival at censored covariate with mean of uncensored covariates before and after.
#'
#' @param at_time A scalar time
#' @param time String column name for the censored covariate.
#' @param delta String column name for the censoring indicator of the covariate.
#' @param surv String column name for the survival estimate.
#' @param data Datafrane containing columns \code{obs}, \code{delta}, and \code{surv}.
#'
#' @return Scalar survival estimate for value \code{at_time}.
#'
#' @export
impute_censored_surv <- function(at_time, time, delta, surv, data) {
  same_time <- which(round(data[, time] - at_time, 8) == 0 & data[, delta] == 1)
  if (length(same_time) == 0) {
    before <- which(data[, time] <= at_time & data[, delta] == 1)
    surv_before <- data[max(before), surv]
    after <- which(data[, time] > at_time & data[, delta] == 1)
    surv_after <- data[min(after), surv]
    return((surv_before + surv_after) / 2)
  } else {
    return(unique(data[same_time, surv][!is.na(data[same_time, surv])]))
  }
}
