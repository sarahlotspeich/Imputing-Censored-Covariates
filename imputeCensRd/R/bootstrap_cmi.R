#' Bootstrap conditional mean imputation (CMI) for a censored predictor
#'
#' Bootstrap standard errors / confidence intervals for conditional mean imputation (CMI) for a censored predictor.
#'
#' @param analysis_model A formula input for the final analysis model of interest. Note that  
#' @param W Column name of observed predictor values (including censored opens). 
#' @param Delta Column name of censoring indicators. Note that \code{Delta = 0} is interpreted as a censored observation. 
#' @param Z Column name of additional fully observed covariates.
#' @param data Dataframe or named matrix containing columns \code{W}, \code{Delta}, and \code{Z}.
#' @param est_surv A string for which CMI approach to be used: fully-parametric (\code{"FP"}), semiparametric (\code{"SP"}), or nonparametric (\code{"NP"}).
#' @param trapezoidal_rule A logical input for whether the trapezoidal rule should be used to approximate the integral in the imputed values. Default is \code{FALSE}.
#' @param dist (If \code{est_surv = "FP"}) The assumed distribution for \code{W} in the AFT model, passed to \code{survival::survreg()}. Default is \code{"weibull"}.
#' @param stratified (If \code{est_surv = "SP"}) If \code{TRUE}, stratification in \code{W} is used to construct time-varying coefficients in the Cox model. Default is \code{FALSE}. 
#' @param surv_between (If \code{est_surv = "SP"} or \code{"NP"}) A string for the method to be used to interpolate for censored values between events. Options include \code{"carry-forward"} (default), \code{"linear"}, or \code{"mean"}.
#' @param surv_beyond (If \code{est_surv = "SP"} or \code{"NP"}) A string for the method to be used to extrapolate the survival curve beyond the last observed event. Options include \code{"drop-off"}, \code{"exponential"} (default), or \code{"weibull"}.
#' @param useSURV (If \code{est_surv = "custom"}) Assumed survival function for \code{W} given \code{Z}. The only arguments to \code{useSURV} should be \code{W} and \code{Z}, in that order.
#' @param B Number of bootstrap resampling replicates. Default is \code{B = 1000}.
#'
#' @return A dataframe containing the following information for the coefficients of \code{analysis_model}:
#' \item{Coeff}{Variable name, corresponding to the user-specified model \code{analysis_model}.}
#' \item{SE}{Bootstrapped standard errors (empirical standard deviation of the bootstrapped estimates).}
#' \item{LB}{Lowerbound of bootstrapped confidence interval (based on quantiles of bootstrapped estimates).}
#' \item{UB}{Upperbound of bootstrapped confidence interval (based on quantiles of bootstrapped estimates).}
#'
#' @export

bootstrap_cmi <- function(analysis_model, W, Delta, Z, data, est_surv, trapezoidal_rule = FALSE, dist = "weibull", stratified = FALSE, surv_between = "carry-forward", surv_beyond = "exponential", useSURV, B = 1000) {
  # Create matrix to hold results from bootstrap replicates 
  re_res <- matrix(data = NA, nrow = B, ncol = (length(Z) + 2))
  
  # Size of resample
  n <- nrow(data)
  
  # Loop through replicates 
  for (b in 1:B) {
    # Sample with replacement from the original data ------------------
    re_rows <- ceiling(runif(n = n, min = 0, max = 1) * nrow(data))
    re_data <- data[re_rows, ]
    
    if (sum(re_data[, Delta]) < n) {
      if (est_surv == "FP") {
        # Use imputeCensRd::cmi_fp() to impute censored x in re_data ------
        re_data_imp <- cmi_fp(W = W, Delta = Delta, Z = Z, data = re_data, 
                              fit = NULL, dist = dist, trapezoidal_rule = trapezoidal_rule)
      } else if (est_surv == "SP") {
        # Use imputeCensRd::cmi_sp() to impute censored x in re_data ------
        re_data_imp <- cmi_sp(W = W, Delta = Delta, Z = Z, data = re_data, 
                              fit = NULL, stratified = stratified, trapezoidal_rule = trapezoidal_rule,
                              surv_between = surv_between, surv_beyond = surv_beyond)
      } else if (est_surv == "NP") {
        # Use imputeCensRd::cmi_np() to impute censored x in re_data ------
        # re_data_imp <- cmi_np(W = W, Delta = Delta, Z = Z, data = re_data, 
        #                      trapezoidal_rule = trapezoidal_rule)
      } else if (est_surv == "custom") {
        # Use imputeCensRd::cmi_custom() to impute censored x in re_data --
        re_data_imp <- cmi_custom(W = W, Delta = Delta, Z = Z, data = re_data, 
                                  useSURV = useSURV, trapezoidal_rule = trapezoidal_rule)
      }
      
      # If imputation was successful, fit the analysis model ------------
      if (re_data_imp$code) {
        re_fit <- lm(formula = analysis_model, data = re_data_imp$imputed_data) 
        ## Save coefficients to results matrix
        re_res[b, ] <- re_fit$coefficients
      }
    } else {
      # If no censored, just fit the usual model
      re_data$imp <- re_data[, W]
      re_fit <- lm(formula = analysis_model, data = re_data) 
      ## Save coefficients to results matrix
      re_res[b, ] <- re_fit$coefficients
    }
  }
  
  # Calculate SE estimate
  se <- apply(X = re_res, MARGIN = 2, FUN = sd, na.rm = TRUE)
  
  # Calculate 95% quantile interval
  lb <- apply(X = re_res, MARGIN = 2, FUN = function(x) quantile(x = x, probs = 0.025, na.rm = TRUE))
  ub <- apply(X = re_res, MARGIN = 2, FUN = function(x) quantile(x = x, probs = 0.975, na.rm = TRUE))
  
  # Return SE and 95% interval
  return(data.frame(Coeff = names(re_fit$coefficients), SE = se, LB = lb, UB = ub))
}