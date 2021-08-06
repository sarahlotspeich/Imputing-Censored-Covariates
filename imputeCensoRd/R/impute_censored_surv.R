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
  # test for bad input
  if (!is.character(time)) { stop("argument time must be a character") }
  if (!is.character(delta)) { stop("argument delta must be a character") }
  if (!is.character(surv)) { stop("argument surv must be a character") }
  if (!is.data.frame(data)) { stop("argument data must be a data frame") }
  # test that data contains columns with specified names
  if (!(time %in% colnames(data))) { stop(paste("data does not have column with name", time)) }
  if (!(delta %in% colnames(data))) { stop(paste("data does not have column with name", delta)) }
  if (!(surv %in% colnames(data))) { stop(paste("data does not have column with name", surv)) }
  # test for improper entries in columns of data
  #### Still deciding whether the following conditions should produce errors (stop()) or warnings
  if (any(data[, time] < 0)) { stop(paste("elements of column", time, "must be positive")) }
  if (!all(data[, delta] %in% c(0, 1))) { stop(paste("elements of column", delta, "must be either 0 or 1")) }
  if (any(data[, surv] < 0 | data[, surv] > 1)) { stop(paste("elements of column", surv, "must be inclusively between 0 and 1"))}
  
  # which (if any) event times are equal to at_time, to 8 decimcal places?
  same_time <- which(round(data[, time] - at_time, 8) == 0 & data[, delta] == 1)
  
  # if no event times are equal to at_time
  if (length(same_time) == 0) {
    # index of greatest event time less than at_time and corresponding survival estimate
    before <- which(data[, time] <= at_time & data[, delta] == 1)
    surv_before <- data[max(before), surv]
   
    # index of smallest event time greater than at_time and corresponding survival estimate
    after <- which(data[, time] > at_time & data[, delta] == 1)
    surv_after <- data[min(after), surv]
    
    # average the above survival estimates
    return((surv_before + surv_after) / 2)
  } else {
    #### SARAH: can you explain what this lien is doing?
    return(unique(data[same_time, surv][!is.na(data[same_time, surv])]))
  }
}
