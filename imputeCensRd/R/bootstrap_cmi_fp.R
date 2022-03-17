#' Bootstrap fully parametric conditional mean imputation (CMI) for a censored predictor
#'
#' Bootstrap fully parametric conditional mean imputation (CMI) for a censored predictor using an accelerated failure-time (AFT) model to estimate conditional survival.
#'
#' @param analysis_model A formula input for the final analysis model of interest. Note that  
#' @param W Column name of observed predictor values (including censored opens). 
#' @param Delta Column name of censoring indicators. Note that \code{Delta = 0} is interpreted as a censored observation. 
#' @param Z Column name of additional fully observed covariates.
#' @param data Dataframe or named matrix containing columns \code{W}, \code{Delta}, and \code{Z}.
#' @param fit A \code{survreg} imputation model object modeling \code{W} on \code{Z}. If \code{fit = NULL} (default), the AFT model with only main effects for \code{Z} and assuming a Weibull distribution is fit internally and used.
#' @param dist (Optional) Assumed distribution for \code{W} in the AFT model, passed to \code{survival::survreg()}. Default is \code{"weibull"}.
#' @param trapezoidal_rule A logical input for whether the trapezoidal rule should be used to approximate the integral in the imputed values. Default is \code{FALSE}.
#' @param B Number of bootstrap resampling replicates. Default is \code{B = 1000}.
#'
#' @return 
#' \item{imputed_data}{A copy of \code{data} with added column \code{imp} containing the imputed values.}
#' \item{code}{Indicator of algorithm status (\code{TRUE} or \code{FALSE}).}
#'
#' @export

bootstrap_cmi_fp <- function(analysis_model, W, Delta, Z, data, fit = NULL, dist = "weibull", trapezoidal_rule = FALSE, B = 1000) {
  # Create matrix to hold results from boostrap replicates 
  re_res <- matrix(data = NA, nrow = B, ncol = (length(Z) + 2))
  
  # Loop through replicates 
  for (b in 1:B) {
    # Sample with replacement from the original data ------------------
    re_rows <- ceiling(runif(n = n, min = 0, max = 1) * nrow(data))
    re_data <- data[re_rows, ]
    
    # Use imputeCensRd::cmi_fp() to impute censored x in re_data ------
    re_data_imp <- cmi_fp(W = W, Delta = Delta, Z = Z, data = re_data, 
                          fit = fit, dist = dist, trapezoidal_rule = trapezoidal_rule)
    
    # If imputation was successful, fit the analysis model ------------
    if (re_data_imp$code) {
      re_fit <- lm(formula = analysis_model, data = re_data_imp$imputed_data) 
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