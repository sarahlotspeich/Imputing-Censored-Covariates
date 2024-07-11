#' Constrained maximum likelihood estimators for the Weibull extensions
#'
#' Constrained maximum likelihood estimators for the Weibull extensions
#'
#' @param imputation_model Imputation model formula (or coercible to formula), a formula expression as for other regression models. The response is usually a survival object as returned by the \code{Surv} function. See the documentation for \code{Surv} for details.
#' @param data Dataframe or named matrix containing columns \code{W}, \code{Delta}, and \code{Z}.
#' @param Xmax (Optional) Upper limit of the domain of the censored predictor. Default is \code{Xmax = Inf}.
#'
#' @return A vector of parameters
#'
#' @export

get_weib_params = function (imputation_model, data, Xmax = Inf) {
  # Extract variable names from imputation_model
  W = all.vars(imputation_model)[1] ## censored covariate
  Delta = all.vars(imputation_model)[2] ## corresponding event indicator
  Z = all.vars(imputation_model)[-c(1:2)] ## additional covariates
  
  # Convert data to dataframe (just in case)
  data = data.frame(data)
  
  # Initialize imputed values 
  data$imp = data[, W] ## start with imp = W
  
  # Fit Cox PH imputation model for X ~ Z 
  fit = coxph(formula = imputation_model, 
              data = data)
  
  # Calculate linear predictor for Cox imputation model or logHRs provided
  lp = predict(fit, reference = "sample") + sum(coef(fit) * fit$means, na.rm = TRUE) ## linear predictors
  
  data$HR = exp(lp)
  be = breslow_estimator(x = NULL, 
                         time = W, 
                         event = Delta, 
                         hr = "HR", 
                         data = data)
  surv_df = with(be, 
                 data.frame(t = times, surv0 = basesurv))

  
  Xtilde = surv_df[nrow(surv_df), "t"]
  SURVmax = surv_df[nrow(surv_df), "surv0"]
  
  if (Xmax < Inf) {
    weibull_params = dbl_constr_weibull(Xtilde = Xtilde, 
                                        rho = SURVmax, 
                                        Xmax = Xmax)
  } else {
    weibull_params = constr_weibull_mle(t = data[, W], 
                                        I_event = data[, Delta], 
                                        Xtilde = Xtilde, 
                                        rho = SURVmax, 
                                        alpha0 = 1e-04)
  }
  
  if (any(is.na(weibull_params))) {
    weibull_params = rep(NA, 2)
  }
  
  return(weibull_params)
}