#' Conditional mean imputation (CMI) for a censored predictor using the Kaplan-Meier estimator
#'
#' Nonparametric conditional mean imputation (CMI) for a censored predictor using a Kaplan-Meier product-limit estimator for estimate survival.
#'
#' @param W Column name of observed predictor values (including censored opens). 
#' @param Delta Column name of censoring indicators. Note that \code{Delta = 0} is interpreted as a censored observation. 
#' @param data Dataframe or named matrix containing columns \code{W}, \code{Delta}, and \code{Z}.
#' @param surv_between A string for the method to be used to interpolate for censored values between events. Options include \code{"carry-forward"} (default), \code{"linear"}, or \code{"mean"}.
#' @param surv_beyond A string for the method to be used to extrapolate the survival curve beyond the last observed event. Options include \code{"drop-off"}, \code{"exponential"} (default), or \code{"weibull"}.
#'
#' @return A copy of \code{data} with added column \code{imp} containing the imputed values.
#'
#' @export
#' @importFrom survival Surv
#' @importFrom survival survfit

cmi_km <- function(W, Delta, data, surv_between = "carry-forward", surv_beyond = "exponential") {
  # Fit the Kaplan-Meier estimator for S(W)
  fit_formula <- as.formula(paste0("Surv(time = ", W, ", event = ", Delta, ") ~ 1"))
  fit <- survfit(formula = fit_formula, 
               data = data)
  surv_df <- data.frame(x = fit$time, surv = fit$surv)
  colnames(surv_df)[1] <- W
  
  # Merge survival estimates into data
  data <- merge(x = data, y = surv_df, all.x = TRUE, sort = FALSE)
  
  # Order data by W
  data <- data[order(data[, W]), ]
  
  # Create an indicator variable for being uncensored
  uncens <- data[, Delta] == 1

  # Make survival at censored W NA so that we can interpolate/extrapolate
  data[!uncens, "surv"] <- NA
  
  # Assume survival at censored W < min(X) = 1
  minX <- data[min(which(uncens)), W] 
  data[which(data[, W] < minX), "surv"] <- 1
  
  # Interpolate baseline survival at censored W < \widetilde{X}
  Xtilde <- data[max(which(uncens)), W] 
  needs_interp <- which(is.na(data[, "surv"]) & data[, W] < Xtilde)
  data[needs_interp, "surv"] <- sapply(X = data[needs_interp, W], 
                                        FUN = interp_surv_between, 
                                        t = surv_df[, W], 
                                        surv = surv_df[, "surv"], 
                                        surv_between = surv_between)
  
  # Extrapolate baseline survival at censored W > \widetilde{X}
  if (surv_beyond == "weibull") {
    # Estimate Weibull parameters using constrained MLE
    SURVmax <- data[max(which(uncens)), "surv"]
    weibull_params <- constr_weibull_mle(t = data[, W], I_event = data[, Delta], Xtilde = Xtilde, rho = SURVmax, alpha0 = 1E-4)
    
    # If weibull params don't converge, quit 
    if (any(is.na(weibull_params))) {
      data$imp <- NA
      return(list(imputed_data = data, code = FALSE))   
    }
    
    # If they do, extrapolate with them 
    needs_extrap <- which(!uncens & data[, W] > Xtilde)
    data[needs_extrap, "surv"] <- sapply(X = data[needs_extrap, W], 
                                          FUN = extrap_surv_beyond, 
                                          t = surv_df[, W], 
                                          surv = surv_df[, "surv"], 
                                          surv_beyond = surv_beyond, 
                                          weibull_params = weibull_params)
  } else {
    needs_extrap <- which(!uncens & data[, W] > Xtilde)
    data[needs_extrap, "surv"] <- sapply(X = data[needs_extrap, W], 
                                          FUN = extrap_surv_beyond, 
                                          t = surv_df[, W], 
                                          surv = surv_df[, "surv"], 
                                          surv_beyond = surv_beyond)
    weibull_params <- NULL
  }
  
  # Calculate imputed values E(X|X>W)
  if (trapezoidal_rule) {
    # Distinct rows (in case of non-unique obs values)
    data_dist <- unique(data[, c(W, Delta, Z, "surv")])
    
    # [T_{(i+1)} - T_{(i)}]
    t_diff <- data_dist[- 1, W] - data_dist[- nrow(data_dist), W]
    
    # S(t+1) + S(t)
    surv_sum <- data_dist[-1, "surv"] + data_dist[- nrow(data_dist), "surv"]
    
    # Use trapezoidal approximation for integral
    for (i in which(!uncens)) {
      sum_surv_i <- sum((data_dist[- nrow(data_dist), W] >= as.numeric(data[i, W])) * surv_sum * t_diff)
      data$imp[i] <- data$imp[i] + (1 / 2) * (sum_surv_i /  data[i, "surv"])
    }
  } else {
    # Builds on the extend_surv function by raising S(t)
    to_integrate <- function(t) {
      surv <- sapply(X = t, 
                     FUN = extend_surv, 
                     t = surv_df[, W], 
                     surv = surv_df[, "surv"], 
                     surv_between = surv_between, 
                     surv_beyond = surv_beyond, 
                     weibull_params = weibull_params)  
      surv
    }
    
    ## Use integrate() to approximate integral from W to \infty of S(t)
    int_surv <- sapply(
      X = which(!uncens), 
      FUN = function(i) { 
        tryCatch(expr = integrate(f = to_integrate, 
                                  lower = data[i, W], 
                                  upper = Inf, 
                                  subdivisions = 2000)$value,
                 error = function(e) return(NA))
      }
    )
    
    ## Calculate E(X|X>W) = W + int_surv / surv(W)
    data$imp[which(!uncens)] <- data[which(!uncens), W] + int_surv / data[which(!uncens), "surv"]
  }
  
  ## Check for infinite/NA imputed values 
  if (any(is.na(data$imp))) {
    data$imp[which(is.na(data$imp))] <- data[which(is.na(data$imp)), W]
  }
  if (any(data$imp  == Inf)) {
    data$imp[which(data$imp == Inf)] <- data[which(data$imp == Inf), W]
  }
  
  # Return input dataset with appended column imp containing imputed values 
  if (any(is.na(data$imp))) {
    return(list(imputed_data = data, code = FALSE))    
  } else {
    return(list(imputed_data = data, code = TRUE))
  }
  
  # Return input dataset with appended column imp containing imputed values 
  return(data)
}