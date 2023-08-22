
#' @export
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
                                   fit = re_imp_mod, 
                                   infinite_integral = infinite_integral)
      
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