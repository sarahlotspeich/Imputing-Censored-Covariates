#' Multiple, semiparametric conditional mean imputation for a censored covariate
#'
#' Multiple, semiparametric conditional mean imputation for a censored covariate using a Cox proportional hazards model and the Breslow's estimator to estimate the conditional survival function.
#'
#' @param imputation_formula imputation model formula (or coercible to formula), a formula expression as for other regression models. The response is usually a survival object as returned by the \code{Surv} function. See the documentation for \code{Surv} for details.
#' @param analysis_formula analysis model formula (or coercible to formula), a formula expression as for other regression models. The response should be a continuous random variable in \code{data}. 
#' @param W Column name of observed (possibly censored) predictor values. 
#' @param Delta Column name of censoring indicators. Note that \code{data[, Delta] = 0} is interpreted as a censored observation. 
#' @param Z Column name of additional fully observed covariates.
#' @param data Dataframe or named matrix containing columns \code{W}, \code{Delta}, and \code{Z}.
#' @param trapezoidal_rule A logical input for whether the trapezoidal rule should be used to approximate the integral in the imputed values. Default is \code{FALSE}.
#' @param Xmax (Optional) Upper limit of the domain of the censored predictor. Default is \code{Xmax = Inf}.
#' @param surv_between A string for the method to be used to interpolate for censored values between events. Options include \code{"cf"} (carry forward, the default), \code{"wm"} (weighted mean), or \code{"m"} (mean).
#' @param surv_beyond A string for the method to be used to extrapolate the survival curve beyond the last observed event. Options include \code{"d"} (immediate drop off), \code{"e"} (exponential extension, the default), or \code{"w"} (weibull extension).
#' @param B numeric, number of imputations. Default is \code{10}. 
#'
#' @return 
#' \item{imputed_data}{A copy of \code{data} with added column \code{imp} containing the imputed values.}
#' \item{code}{Indicator of algorithm status (\code{TRUE} or \code{FALSE}).}
#'
#' @export

cmi_sp_bootstrap = function(imputation_formula, analysis_formula, W, Delta, Z, data, trapezoidal_rule = FALSE, Xmax = Inf, surv_between = "cf", surv_beyond = "e", maxiter = 100, B = 10) {
  # Size of resample
  n = nrow(data)
  
  # Loop through replicates 
  for (b in 1:B) {
    # Sample with replacement from the original data ------------------
    re_rows = ceiling(runif(n = n, min = 0, max = 1) * nrow(data))
    re_data = data[re_rows, ]
    
    # Check for censoring in resampled data
    if (sum(re_data[, Delta]) < n) {
      # Fit AFT imputation model for X ~ Z 
      re_imp_mod = survreg(formula = imputation_formula, 
                           data = re_data, 
                           dist = "weibull")
      
      # Use imputeCensRd::cmi_fp_weibull() to impute censored x in re_data ------
      re_data_imp = cmi_sp(W = W, 
                           Delta = Delta, 
                           Z = Z, 
                           data = re_data, 
                           fit = re_imp_mod, 
                           trapezoidal_rule = trapezoidal_rule, 
                           Xmax = Xmax,
                           surv_between = surv_between, 
                           surv_beyond = surv_beyond)
      
      # If imputation was successful, fit the analysis model ------------
      if (re_data_imp$code) {
        re_fit = lm(formula = analysis_formula, 
                    data = re_data_imp$imputed_data)
      }
    } else {
      # If no censored, just fit the usual model
      re_data$imp = re_data[, W]
      re_fit = lm(formula = analysis_model, data = re_data) 
    }
    
    ## Create matrix to hold results from bootstrap replicates 
    if (b == 1) {
      re_res = re_var = matrix(data = NA, 
                               nrow = B, 
                               ncol = length(re_fit$coefficients))
    }
    
    ## Save coefficients to results matrix
    re_res[b, ] = re_fit$coefficients
    re_var[b, ] = diag(vcov(re_fit))
  } 
  
  # Return table of pooled estimates
  tab = data.frame(Coefficient = names(re_fit$coefficients),
                   Est = colMeans(re_res),
                   SE = sqrt(colMeans(re_var) + (B + 1) * colMeans((re_res - colMeans(re_res)) ^ 2)))
  return(tab)  
}