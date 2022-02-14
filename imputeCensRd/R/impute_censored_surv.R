impute_censored_surv <- function(at_time, time, event, surv, data) {
  # which (if any) event times are equal to at_time?
  same_time <- which(round(data[, time] - at_time, 8) == 0 & data[, event] == 1 & !is.na(data[, surv]))
  
  # if no event times are equal to at_time, impute with the mean of values immediately before/after
  if (length(same_time) == 0) {
    # index of greatest event time less than at_time
    before <- which(data[, time] <= at_time & data[, event] == 1)
    # corresponding survival estimate
    surv_before <- data[max(before), surv]
    
    # index of smallest event time greater than at_time
    after <- which(data[, time] > at_time & data[, event] == 1)
    # corresponding survival estimate
    surv_after <- data[min(after), surv]
    
    # linear interpolation of survival estimates before and after
    return(surv_before + (surv_after - surv_before) / (data[min(after), time] - data[max(before), time]) * (at_time - data[max(before), time]))
  } else {
    return(data[max(same_time), surv])
  }
}