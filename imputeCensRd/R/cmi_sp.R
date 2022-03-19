#' Semiparametric conditional mean imputation (CMI) for a censored predictor
#'
#' Semiparametric conditional mean imputation (CMI) for a censored predictor using a Cox proportional hazards model and the Breslow estimator to estimate conditional survival.
#'
#' @param W Column name of observed predictor values (including censored opens). 
#' @param Delta Column name of censoring indicators. Note that \code{data[, Delta] = 0} is interpreted as a censored observation. 
#' @param Z Column name of additional fully observed covariates.
#' @param data Dataframe or named matrix containing columns \code{W}, \code{Delta}, and \code{Z}.
#' @param fit A \code{coxph} imputation model object modeling \code{W} on \code{Z}. If \code{fit = NULL} (default), the Cox model with only main effects for \code{Z} is fit internally and used.
#' @param trapezoidal_rule A logical input for whether the trapezoidal rule should be used to approximate the integral in the imputed values. Default is \code{FALSE}.
#' @param surv_between A string for the method to be used to interpolate for censored values between events. Options include \code{"carry-forward"} (default), \code{"linear"}, or \code{"mean"}.
#' @param surv_beyond A string for the method to be used to extrapolate the survival curve beyond the last observed event. Options include \code{"carry-forward"} (default), \code{"drop-off"}, \code{"exponential"}, or \code{"weibull"}.
#' @param force_last_event A logical input to force the last observed value for \code{W} to be treated as an event. Default is \code{FALSE}.
#'
#' @return 
#' \item{imputed_data}{A copy of \code{data} with added column \code{imp} containing the imputed values.}
#' \item{code}{Indicator of algorithm status (\code{TRUE} or \code{FALSE}).}
#'
#' @export
#' @importFrom survival coxph Surv

cmi_sp <- function(W, Delta, Z, data, fit = NULL, trapezoidal_rule = FALSE, surv_between = "carry-forward", surv_beyond = "carry-forward", force_last_event = FALSE) {
  # Assume last observed value is an event regardless
  if (force_last_event) {
    data[which.max(data[, W]), Delta] <- 1
  }
  
  # If no imputation model was supplied, fit a Cox PH using main effects
  if (is.null(fit)) {
    fit_formula <- as.formula(paste0("Surv(time = ", W, ", event = ", Delta, ") ~ ", paste0(Z, collapse = " + ")))
    fit <- coxph(formula = fit_formula, 
                 data = data)
  }
  
  # Calculate linear predictor \lambda %*% Z for Cox model
  lp <- data.matrix(data[, Z]) %*% matrix(data = fit$coefficients, ncol = 1)
  
  ## Calculate hazard ratio for Cox model
  data$HR <- exp(lp)
  
  # Estimate baseline survival from Cox model fit using Breslow's estimator
  be <- breslow_estimator(x = NULL, 
                          time = W, 
                          event = Delta, 
                          hr = "HR", 
                          data = data)
  surv_df <- with(be, data.frame(t = times, surv0 = basesurv))
  colnames(surv_df)[1] <- W
  ## Merge baseline survival estimates into data
  data <- merge(x = data, y = surv_df, all.x = TRUE, sort = FALSE)
  
  # Order data by W
  data <- data[order(data[, W]), ]
  
  # Create an indicator variable for being uncensored
  uncens <- data[, Delta] == 1
  
  # For people with events, obs = X
  data$imp <- data[, W]
  
  # Assume survival at censored W < min(X) = 1
  minX <- data[min(which(uncens)), W] 
  data[which(data[, W] < minX), "surv0"] <- 1
  
  # Interpolate baseline survival at censored W < \widetilde{X}
  Xtilde <- data[max(which(uncens)), W] 
  needs_interp <- which(is.na(data[, "surv0"]) & data[, W] < Xtilde)
  data[needs_interp, "surv0"] <- sapply(X = data[needs_interp, W], 
                                        FUN = interp_surv_between, 
                                        t = surv_df[, W], 
                                        surv = surv_df[, "surv0"], 
                                        surv_between = surv_between)
  
  # Extrapolate baseline survival at censored W > \widetilde{X}
  if (surv_beyond == "weibull") {
    # Estimate Weibull parameters using constrained MLE
    SURVmax <- data[max(which(uncens)), "surv0"]
    weibull_params <- constr_weibull_mle(t = data[, W], I_event = data[, Delta], Xtilde = Xtilde, rho = SURVmax, alpha0 = 0.1)
    
    # If weibull params don't converge, quit 
    if (any(is.na(weibull_params))) {
      data$imp <- NA
      return(list(imputed_data = data, code = FALSE))   
    }
    
    # If they do, extrapolate with them 
    needs_extrap <- which(!uncens & data[, W] > Xtilde)
    data[needs_extrap, "surv0"] <- sapply(X = data[needs_extrap, W], FUN = extrap_surv_beyond, 
                                          t = surv_df[, W], surv = surv_df[, "surv0"], surv_beyond = surv_beyond, weibull_params = weibull_params)
  } else {
    needs_extrap <- which(!uncens & data[, W] > Xtilde)
    data[needs_extrap, "surv0"] <- sapply(X = data[needs_extrap, W], 
                                          FUN = extrap_surv_beyond, 
                                          t = surv_df[, W], 
                                          surv = surv_df[, "surv0"], 
                                          surv_beyond = surv_beyond)
    weibull_params <- NULL
  }

  # Calculate conditional survival W | Z for all 
  data$surv <- data[, "surv0"] ^ data[, "HR"]
  
  # Calculate imputed values E(X|X>W,Z)
  if (trapezoidal_rule) {
    # Distinct rows (in case of non-unique obs values)
    data_dist <- unique(data[, c(W, Delta, Z, "surv")])
    
    # [T_{(i+1)} - T_{(i)}]
    t_diff <- data_dist[- 1, W] - data_dist[- nrow(data_dist), W]
    
    # S(t+1|Z) + S(t|Z)
    surv_sum <- data_dist[-1, "surv"] + data_dist[- nrow(data_dist), "surv"]
    
    # Censored subject values W (to be imputed)
    t_cens <- data[which(!uncens), W]
    
    # Use trapezoidal approximation for integral
    for (i in which(!uncens)) {
      sum_surv_i <- sum((data_dist[-nrow(data_dist), W] >= as.numeric(data[i, "surv"])) * surv_sum * t_diff)
      data$imp[i] <- data$imp[i] + (1 / 2) * (sum_surv_i /  data[i, "surv"])
    }
  } else {
    # Builds on the extend_surv function by raising S_0(t) ^ HR = S(t|Z)
    to_integrate <- function(t, hr) {
      basesurv <- sapply(X = t, 
                         FUN = extend_surv, 
                         t = surv_df[, W], 
                         surv = surv_df[, "surv0"], 
                         surv_between = surv_between, 
                         surv_beyond = surv_beyond, 
                         weibull_params = weibull_params)  
      basesurv ^ as.numeric(hr)
    }
    
    ## Use integrate() to approximate integral from W to \infty of S(t|Z)
    int_surv <- sapply(
      X = which(!uncens), 
      FUN = function(i) { 
        tryCatch(expr = integrate(f = to_integrate, 
                                  lower = data[i, W], 
                                  upper = Inf, 
                                  subdivisions = 2000, 
                                  hr = data[i, "HR"])$value,
                 error = function(e) return(NA))
      }
    )
    
    ## Calculate E(X|X>W,Z) = W + int_surv / surv(W|Z)
    data$imp[which(!uncens)] <- data[which(!uncens), W] + int_surv / data[which(!uncens), "surv"]
  }
  
  ## Check for infinite imputed values 
  if (any(data$imp  == Inf)) {
    data$imp[which(data$imp == Inf)] <- data[which(data$imp == Inf), W]
  }
  
  # Return input dataset with appended column imp containing imputed values 
  if (any(is.na(data$imp))) {
    return(list(imputed_data = data, code = FALSE))  
  } else {
    return(list(imputed_data = data, code = TRUE))
  }
}