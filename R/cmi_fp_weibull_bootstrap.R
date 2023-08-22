#' Multiple, fully parametric conditional mean imputation for a censored covariate (Weibull distribution)
#'
#' Multiple, fully parametric conditional mean imputation for a censored covariate using an accelerated failure-time model with a Weibull distribution to estimate the conditional survival function.
#'
#' @param imputation_formula imputation model formula (or coercible to formula), a formula expression as for other regression models. The response is usually a survival object as returned by the \code{Surv} function. See the documentation for \code{Surv} for details.
#' @param analysis_formula analysis model formula (or coercible to formula), a formula expression as for other regression models. The response should be a continuous random variable in \code{data}. 
#' @param W character, column name for observed values of the censored covariate 
#' @param Delta character, column name for censoring indicators. Note that \code{Delta = 0} is interpreted as a censored observation. 
#' @param Z character vector, column name(s) of additional fully observed covariate(s).
#' @param data Dataframe or named matrix containing columns \code{W}, \code{Delta}, and \code{Z}.
#' @param infinite_integral (optional) logical, if \code{infinite_integral = TRUE} (default) then conditional means are found by integrating from \code{W} to \code{Inf}, whereas if \code{infinite_integral = FALSE} they are found by subtracting the integral from \code{0} to \code{W} from the mean. If \code{infinite_integral = NA} instead, the analytical solutions are used to find the conditional means.
#' @param max_iter (optional) numeric, maximum iterations allowed in call to \code{survival::survreg()}. Default is \code{100}.
#' @param B numeric, number of imputations. Default is \code{10}. 
#'
#' @return 
#' \item{imputed_data}{A copy of \code{data} with added column \code{imp} containing the imputed values.}
#' \item{code}{Indicator of algorithm status (\code{TRUE} or \code{FALSE}).}
#'
#' @export
#' @importFrom survival survreg 
#' @importFrom survival Surv 
#' @import survival

cmi_fp_weibull_bootstrap = function(imputation_formula, analysis_formula, W, Delta, Z, data, infinite_integral = TRUE, maxiter = 100, B = 10) {
  # Size of resample
  n = nrow(data)
  
  # Loop through replicates 
  for (b in 1:B) {
    # Sample with replacement from the original data ------------------
    re_rows = ceiling(runif(n = n, min = 0, max = 1) * nrow(data))
    re_data = data[re_rows, ]
    
    # Check for censoring in resampled data
    if (sum(re_data[, Delta]) < n) {
      # Use imputeCensRd::cmi_fp_weibull() to impute censored x in re_data ------
      re_data_imp = cmi_fp_weibull(imputation_formula = imputation_formula, 
                                   W = W, 
                                   Delta = Delta, 
                                   Z = Z, 
                                   data = re_data, 
                                   infinite_integral = infinite_integral)
      
      # If imputation was successful, fit the analysis model ------------
      if (re_data_imp$code) {
        re_fit = lm(formula = analysis_formula, 
                    data = re_data_imp$imputed_data)
      }
    } else {
      # If no censored, just fit the usual model
      re_data$imp = re_data[, W]
      re_fit = lm(formula = analysis_model, 
                  data = re_data) 
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