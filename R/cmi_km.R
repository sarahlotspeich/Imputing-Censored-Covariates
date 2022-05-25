#' Conditional mean imputation (CMI) for a censored predictor using the Kaplan-Meier estimator
#'
#' Nonparametric conditional mean imputation (CMI) for a censored predictor using a Kaplan-Meier product-limit estimator for estimate survival.
#'
#' @param W Column name of observed predictor values (including censored opens). 
#' @param Delta Column name of censoring indicators. Note that \code{Delta = 0} is interpreted as a censored observation. 
#' @param Z Column name of additional fully observed categorical covariate. If provided and \code{stratified = TRUE}, the Kaplan-Meier estimator will be stratified on \code{Z}.
#' @param data Dataframe or named matrix containing columns \code{W}, \code{Delta}, and \code{Z}.
#' @param stratified A logical input for whether the Kaplan-Meier estimator should be stratified on \code{Z}. Default is \code{TRUE}.
#' @param trapezoidal_rule A logical input for whether the trapezoidal rule should be used to approximate the integral in the imputed values. Default is \code{FALSE}.
#' @param Xmax (Optional) Upper limit of the domain of the censored predictor. Default is \code{Xmax = Inf}.
#' @param surv_between A string for the method to be used to interpolate for censored values between events. Options include \code{"carry-forward"} (default), \code{"linear"}, or \code{"mean"}.
#' @param surv_beyond A string for the method to be used to extrapolate the survival curve beyond the last observed event. Options include \code{"drop-off"}, \code{"exponential"} (default), or \code{"weibull"}.
#'
#' @return A copy of \code{data} with added column \code{imp} containing the imputed values.
#'
#' @export
#' @importFrom survival Surv
#' @importFrom survival survfit

cmi_km <- function(W, Delta, Z = NULL, data, stratified = TRUE, trapezoidal_rule = FALSE, Xmax = Inf, surv_between = "carry-forward", surv_beyond = "exponential") {
  # Fit the Kaplan-Meier estimator for S(W)
  fit_formula <- as.formula(paste0("Surv(time = ", W, ", event = ", Delta, ") ~ 1"))
  
  # If (categorical) Z is supplied, fit separate models for Z = 0/ Z = 1
  if (!is.null(Z) & stratified) {
    levelsZ <- unique(data[, Z])
    surv_df <- data.frame()
    
    for (z in levelsZ) {
      # Fit KM within strata with Z = z
      has_z <- data[, Z] == z
      fit_z <- survfit(formula = fit_formula, 
                     data = data[has_z, ])
      surv_df_z <- data.frame(x = fit_z$time, z, surv = fit_z$surv)
      colnames(surv_df_z)[1:2] <- c(W, Z)
      
      # Combine 
      surv_df <- rbind(surv_df, surv_df_z)
    }
  } else {
    fit <- survfit(formula = fit_formula, 
                   data = data)
    surv_df <- data.frame(x = fit$time, surv = fit$surv)
    colnames(surv_df)[1] <- W
  }
  
  # Merge survival estimates into data
  data <- merge(x = data, 
                y = surv_df, 
                all.x = TRUE, 
                sort = FALSE)
  
  # Order data by W
  data <- data[order(data[, W]), ]
  
  # Save data frame with KM estimates (ordered)
  km_surv <- data[data[, Delta] == 1, ]
  
  # Create an indicator variable for being uncensored
  uncens <- data[, Delta] == 1

  # Make survival at censored W NA so that we can interpolate/extrapolate
  data[!uncens, "surv"] <- NA
  
  # For people with events, obs = X
  data$imp <- data[, W]
  
  # If stratifying on Z, do the following within strata Z = 0 / Z = 1
  if (!is.null(Z) & stratified) {
    for (z in levelsZ) {
      # Indicator for having Z = z
      has_z <- data[, Z] == z
      
      # Subset to data having Z = z
      data_z <- data[has_z, ]
      
      # Extend survival curve and impute censored W
      res <- extend_and_impute(data = data_z, 
                               W = W, 
                               Delta = Delta, 
                               Z = Z, 
                               S = "surv", 
                               S0 = NULL, 
                               Xmax = Xmax, 
                               trapezoidal_rule = trapezoidal_rule, 
                               surv_between = surv_between, 
                               surv_beyond = surv_beyond)
      
      # Combine imputed data_z with the rest of data 
      if (res$code) {
        data <- rbind(data[!has_z, ], res$data)
      } else {
        data$imp <- NA
        return(list(imputed_data = data, code = FALSE))
      }
    }
  } else {
    # Extend survival curve and impute censored W
    res <- extend_and_impute(data = data, 
                             W = W, 
                             Delta = Delta, 
                             Z = Z, 
                             S = "surv", 
                             S0 = NULL, 
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
  }
  
  ## Check for infinite/NA imputed values 
  if (any(is.na(data$imp))) {
    data$imp[which(is.na(data$imp))] <- data[which(is.na(data$imp)), W]
  }
  if (any(data$imp  == Inf)) {
    data$imp[which(data$imp == Inf)] <- data[which(data$imp == Inf), W]
  }
  
  # Return input dataset with appended column imp containing imputed values 
  return(list(imputed_data = data, code = !any(is.na(data$imp))))
}
