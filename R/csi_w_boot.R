#' Semiparametric single conditional mean imputation (CMI) with boostrapping for a censored predictor with
#'
#' Semiparametric single conditional mean imputation (CMI) with bootstrapping for a censored predictor using a Cox proportional hazards model and the Breslow estimator to estimate conditional survival.
#'
#' @param imputation_model Imputation model formula (or coercible to formula), a formula expression as for other regression models. The response is usually a survival object as returned by the \code{Surv} function. See the documentation for \code{Surv} for details.
#' @param analysis_model Analysis model formula (or coercible to formula), a formula expression as for other regression models. The response should be a continuous outcome for normal linear regression.
#' @param data Dataframe or named matrix containing columns \code{W}, \code{Delta}, and \code{Z}.
#' @param trapezoidal_rule A logical input for whether the trapezoidal rule should be used to approximate the integral in the imputed values. Default is \code{FALSE}.
#' @param Xmax (Optional) Upper limit of the domain of the censored predictor. Default is \code{Xmax = Inf}.
#' @param surv_between A string for the method to be used to interpolate for censored values between events. Options include \code{"cf"} (carry forward, the default), \code{"wm"} (weighted mean), or \code{"m"} (mean).
#' @param surv_beyond A string for the method to be used to extrapolate the survival curve beyond the last observed event. Options include \code{"d"} (immediate drop off), \code{"e"} (exponential extension, the default), or \code{"w"} (weibull extension).
#' @param B numeric, number of bootstraps used for standard errors Default is \code{500}. 
#'
#' @return a \code{dataframe} of pooled coefficient and standard error estimates 
#'
#' @export

csi_w_boot = function (imputation_model, analysis_model, data, trapezoidal_rule = FALSE, Xmax = Inf, surv_between = "cf", surv_beyond = "e", B = 500) {
  # Impute censored covariates in the original data 
  orig_imp = cmi_sp(imputation_model = imputation_model, 
                    lp = NULL, 
                    data = data, 
                    trapezoidal_rule = trapezoidal_rule, 
                    Xmax = Xmax, 
                    surv_between = surv_between, 
                    surv_beyond = surv_beyond)
  
  # Extract variable names from imputation_model
  W = all.vars(imputation_model)[1] ## censored covariate
  
  # Take names and dimension from naive fit 
  data[, "imp"] = data[, W] ## start imputed value with observed 
  naive_fit = lm(formula = analysis_model, 
                 data = data) 
  p = length(naive_fit$coefficients) ## dimension of beta vector 

  # Initialize empty matrix to hold coefficient and variance estimates 
  beta_b = vbeta_b = matrix(data = NA, 
                            nrow = B, 
                            ncol = p)
  
  # Check for successful imputation 
  if (orig_imp$code) {
    ## Bootstrap
    for (b in 1:B) {
      ### Resample with replacement from orig_imp$imputed_data -----------------
      re_rows = ceiling(runif(n = n, min = 0, max = 1) * nrow(data))
      re_data = orig_imp$imputed_data[re_rows, ]
      
      ### Fit analysis_model to resampled data ---------------------------------
      re_fit = lm(formula = analysis_model, 
                  data = re_data) 
      
      ### Save coefficient and variances estimates -----------------------------
      beta_b[b, ] = coefficients(re_fit)
      vbeta_b[b, ] = diag(vcov(re_fit))
    }
    
    ## Pool coefficient and variance estimates 
    beta_pooled = colMeans(beta_b) 
    beta_pooled_rep = matrix(data = beta_pooled, 
                             nrow = B, 
                             ncol = p, 
                             byrow = TRUE)
    vbeta_within = colMeans(vbeta_b)
    vbeta_between = colSums((beta_b - beta_pooled_rep) ^ 2) / (B - 1) 
    vbeta_pooled = vbeta_within + vbeta_between 
    
    ## Construct table of results 
    tab = data.frame(Coefficient = names(naive_fit$coefficients),
                     Est = beta_pooled,
                     SE = sqrt(vbeta_pooled))
  } else {
    ## If imputation unsuccessful, return empty table of NA results
    tab = data.frame(Coefficient = names(naive_fit$coefficients),
                     Est = NA,
                     SE = NA)
  }
  
  # Return table of pooled estimates 
  return(tab)  
}