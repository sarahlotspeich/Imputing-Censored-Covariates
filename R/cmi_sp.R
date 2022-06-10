#' Semiparametric conditional mean imputation (CMI) for a censored predictor
#'
#' Semiparametric conditional mean imputation (CMI) for a censored predictor using a Cox proportional hazards model and the Breslow estimator to estimate conditional survival.
#'
#' @param W Column name of observed (possibly censored) predictor values. 
#' @param Delta Column name of censoring indicators. Note that \code{data[, Delta] = 0} is interpreted as a censored observation. 
#' @param Z Column name of additional fully observed covariates.
#' @param data Dataframe or named matrix containing columns \code{W}, \code{Delta}, and \code{Z}.
#' @param trapezoidal_rule A logical input for whether the trapezoidal rule should be used to approximate the integral in the imputed values. Default is \code{FALSE}.
#' @param Xmax (Optional) Upper limit of the domain of the censored predictor. Default is \code{Xmax = Inf}.
#' @param surv_between A string for the method to be used to interpolate for censored values between events. Options include \code{"carry-forward"} (default), \code{"linear"}, or \code{"mean"}.
#' @param surv_beyond A string for the method to be used to extrapolate the survival curve beyond the last observed event. Options include \code{"drop-off"}, \code{"exponential"} (default), or \code{"weibull"}.
#'
#' @return 
#' \item{imputed_data}{A copy of \code{data} with added column \code{imp} containing the imputed values.}
#' \item{code}{Indicator of algorithm status (\code{TRUE} or \code{FALSE}).}
#'
#' @export
#' @importFrom survival coxph 
#' @importFrom survival Surv

cmi_sp <- function(W, Delta, Z, data, trapezoidal_rule = FALSE, Xmax = Inf, surv_between = "carry-forward", surv_beyond = "exponential") {
  # Fit a Cox PH using main effects
  fit_formula <- as.formula(paste0("Surv(time = ", W, ", event = ", Delta, ") ~ ", paste0(Z, collapse = " + ")))
  fit <- coxph(formula = fit_formula, 
               data = data)
  
  # Calculate linear predictor \lambda %*% Z for Cox model
  lp <- predict(fit, reference="sample") + sum(coef(fit) * fit$means, na.rm = TRUE)
  
  ## Calculate hazard ratio for Cox model
  data$HR <- exp(lp)
  
  # Estimate baseline survival from Cox model fit using Breslow's estimator
  be <- breslow_estimator(x = NULL, 
                          time = W, 
                          event = Delta, 
                          hr = "HR", 
                          data = data)
  surv_df <- with(be, 
                  data.frame(t = times, 
                             surv0 = basesurv))
  colnames(surv_df)[1] <- W
  
  ## Merge baseline survival estimates into data
  data <- merge(x = data, 
                y = surv_df, 
                all.x = TRUE, 
                sort = FALSE)
  
  # Order data by W
  data <- data[order(data[, W]), ]
  
  # Create an indicator variable for being uncensored
  uncens <- data[, Delta] == 1
  
  # For people with events, obs = X
  data$imp <- data[, W]
  
  # Extend survival curve and impute censored W
  res <- extend_and_impute_cox(data = data, 
                               W = W, 
                               Delta = Delta, 
                               Z = Z, 
                               S0 = "surv0", 
                               Xmax = Xmax, 
                               trapezoidal_rule = trapezoidal_rule, 
                               surv_between = surv_between, 
                               surv_beyond = surv_beyond)
  if (res$code) {
    data <- res$data
  } else {
    data$imp <- NA
    return(list(imputed_data = data, code = FALSE))
  }

  ## Check for infinite/NA imputed values 
  if (any(is.na(data$imp))) {
    data$imp[which(is.na(data$imp))] <- data[which(is.na(data$imp)), W]
  }
  if (any(data$imp  == Inf)) {
    data$imp[which(data$imp == Inf)] <- data[which(data$imp == Inf), W]
  }
  
  # Return input dataset with appended column imp containing imputed values 
  return(list(imputed_data = data, 
              code = !any(is.na(data$imp)))) 
}