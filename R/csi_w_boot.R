#' Semiparametric single conditional mean imputation (CMI) with boostrapping for a censored predictor with
#'
#' Semiparametric single conditional mean imputation (CMI) with bootstrapping for a censored predictor using a Cox proportional hazards model and the Breslow estimator to estimate conditional survival.
#'
#' @param imputation_model Imputation model formula (or coercible to formula), a formula expression as for other regression models. The response is usually a survival object as returned by the \code{Surv} function. See the documentation for \code{Surv} for details.
#' @param analysis_model Analysis model formula (or coercible to formula), a formula expression as for other regression models. The response should be a continuous outcome for normal linear regression.
#' @param data Dataframe or named matrix containing columns \code{W}, \code{Delta}, and \code{Z}.
#' @param integral A string input for how to approximate the integral in the imputed values. Default is \code{integral="aq"} for adaptive quadrature, but \code{"tr"} (trapezoidal rule) is also available.
#' @param Xmax (Optional) Upper limit of the domain of the censored predictor. Default is \code{Xmax = Inf}.
#' @param subdivisions (Optional) Passed through to \code{integrate}, the maximum number of subintervals. Default is \code{subdivisions = 2000}.
#' @param surv_between A string for the method to be used to interpolate for censored values between events. Options include \code{"cf"} (carry forward, the default), \code{"wm"} (weighted mean), or \code{"m"} (mean).
#' @param surv_beyond A string for the method to be used to extrapolate the survival curve beyond the last observed event. Options include \code{"d"} (immediate drop off), \code{"e"} (exponential extension, the default), or \code{"w"} (weibull extension).
#' @param B Numeric, number of bootstraps used for standard errors Default is \code{500}. 
#' @param stratify Logical, if \code{stratify = FALSE} (the default), the bootstraps resample from the pooled sample. If \code{stratify = TRUE}, the bootstraps resample separately from the censored and uncensored subsamples.
#'
#' @return a \code{dataframe} of pooled coefficient and standard error estimates 
#'
#' @export

csi_w_boot = function (imputation_model, analysis_model, data, integral = "aq", Xmax = Inf, subdivisions = 2000, surv_between = "cf", surv_beyond = "e", B = 500, stratify = FALSE) {
  # Impute censored covariates in the original data 
  orig_imp = cmi_sp(imputation_model = imputation_model, 
                    lp = NULL, 
                    data = data, 
                    integral = integral,
                    Xmax = Xmax, 
                    subdivisions = subdivisions, 
                    surv_between = surv_between, 
                    surv_beyond = surv_beyond)
  
  # Extract variable names from imputation_model
  W = all.vars(imputation_model)[1] ## censored covariate
  Delta = all.vars(imputation_model)[2] ## corresponding event indicator
  
  # Define separate objects for censored and uncensored observations -----------
  uncens_data = orig_imp$imputed_data[orig_imp$imputed_data[, Delta] == 1, ]
  cens_data = orig_imp$imputed_data[orig_imp$imputed_data[, Delta] == 0, ]
  
  # Check for successful imputation 
  if (orig_imp$code) {
    ## Take names and dimension from single imputation fit 
    csi_fit = lm(formula = analysis_model, 
                 data = orig_imp$imputed_data) 
    p = length(csi_fit$coefficients) ## dimension of beta vector 
    
    ## Initialize empty matrix to hold coefficient and variance estimates 
    beta_b = vbeta_b = matrix(data = NA, 
                              nrow = B, 
                              ncol = p)
    ## Bootstrap
    for (b in 1:B) {
      ## Resample with replacement from orig_imp$imputed_data -----------------
      if (stratify) {
        ## Uncensored rows ---------------------------------------------------------
        re_uncens_rows = ceiling(runif(n = nrow(uncens_data), min = 0, max = 1) * nrow(uncens_data))
        re_uncens_data = uncens_data[re_uncens_rows, ]
        ## Censored rows -----------------------------------------------------------
        re_cens_rows = ceiling(runif(n = nrow(cens_data), min = 0, max = 1) * nrow(cens_data))
        re_cens_data = cens_data[re_cens_rows, ]
        ## Put them together -------------------------------------------------------
        re_data = rbind(re_uncens_data, 
                        re_cens_data)
      } else {
        re_rows = ceiling(runif(n = nrow(data), min = 0, max = 1) * nrow(data))
        re_data = orig_imp$imputed_data[re_rows, ]
      }
      
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
    tab = data.frame(Coefficient = names(csi_fit$coefficients),
                     NEst = csi_fit$coefficients, 
                     BEst = beta_pooled,
                     NSE = sqrt(diag(vcov(csi_fit))),
                     BSE = sqrt(vbeta_pooled))
  } else {
    ## Take names from naive fit 
    data[, "imp"] = data[, W]
    naive_fit = lm(formula = analysis_model, 
                   data = data) 
    
    ## If imputation unsuccessful, return empty table of NA results
    tab = data.frame(Coefficient = names(naive_fit$coefficients),
                     NEst = NA,
                     BEst = NA,
                     NSE = NA, 
                     BSE = NA)
  }
  
  # Return table of pooled estimates 
  return(tab)  
}