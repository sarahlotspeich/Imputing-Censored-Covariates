#' Semiparametric conditional mean imputation (CMI) for a censored predictor
#'
#' Semiparametric conditional mean imputation (CMI) for a censored predictor using a Cox proportional hazards model and the Breslow estimator to estimate conditional survival.
#'
#' @param imputation_model Imputation model formula (or coercible to formula), a formula expression as for other regression models. The response is usually a survival object as returned by the \code{Surv} function. See the documentation for \code{Surv} for details.
#' @param data Dataframe or named matrix containing columns \code{W}, \code{Delta}, and \code{Z}.
#' @param trapezoidal_rule A logical input for how to approximate the integral in the imputed values. Default is \code{trapezoidal_rule=FALSE} to use adaptive quadrature, but \code{"tr"}, and \code{trapezoidal_rule=TRUE} will use the trapezoidal rule.
#' @param Xmax (Optional) Upper limit of the domain of the censored predictor. Default is \code{Xmax = Inf}.
#' @param subdivisions (Optional) Passed through to \code{integrate}, the maximum number of subintervals. Default is \code{subdivisions = 2000}.
#' @param surv_between A string for the method to be used to interpolate for censored values between events. Options include \code{"cf"} (carry forward, the default), \code{"wm"} (weighted mean), or \code{"m"} (mean).
#' @param surv_beyond A string for the method to be used to extrapolate the survival curve beyond the last observed event. Options include \code{"d"} (immediate drop off), \code{"e"} (exponential extension, the default), or \code{"w"} (weibull extension).
#'
#' @return 
#' \item{imputed_data}{A copy of \code{data} with added column \code{imp} containing the imputed values.}
#' \item{code}{Indicator of algorithm status (\code{TRUE} or \code{FALSE}).}
#'
#' @export

cmi_sp = function (imputation_model, data, trapezoidal_rule = TRUE, Xmax = Inf, subdivisions = 2000, surv_between = "cf", surv_beyond = "e") {
  # Extract variable names from imputation_model
  W = all.vars(imputation_model)[1] ## censored covariate
  Delta = all.vars(imputation_model)[2] ## corresponding event indicator
  Z = all.vars(imputation_model)[-c(1:2)] ## additional covariates
  
  # Convert data to dataframe (just in case)
  data = data.frame(data)
  
  # Initialize imputed values 
  data$imp = data[, W] ## start with imp = W
  
  # Fit Cox PH imputation model for X ~ Z 
  fit = tryCatch(expr = {
    coxph(formula = imputation_model, 
          data = data)}, 
    warning = function(w) {
      list(coefficients = rep(NA, length(Z)))
    }
  )
  
  # Check for NA values in logHRs 
  if (any(is.na(fit$coefficients))) {
    data$imp = NA ## make imp = NA
    return(list(imputed_data = data, 
                code = !any(is.na(data$imp))))
  }
  
  # Calculate linear predictor for Cox imputation model
  lp = predict(fit, reference = "sample") + sum(coef(fit) * fit$means, na.rm = TRUE) ## linear predictors
  
  data$HR = exp(lp)
  be = breslow_estimator(x = NULL, 
                         time = W, 
                         event = Delta, 
                         hr = "HR", 
                         data = data)
  surv_df = with(be, 
                  data.frame(t = times, surv0 = basesurv))
  colnames(surv_df)[1] = W
  data = merge(x = data, 
                y = surv_df, 
                all.x = TRUE, 
                sort = FALSE)
  data = data[order(data[, W]), ]
  uncens = data[, Delta] == 1
  minX = data[min(which(uncens)), W]
  data[which(data[, W] < minX), "surv0"] = 1
  Xtilde = data[max(which(uncens)), W]
  needs_interp = which(is.na(data[, "surv0"]) & data[, W] < Xtilde)
  data[needs_interp, "surv0"] = sapply(X = data[needs_interp, W], 
                                        FUN = interp_surv_between, 
                                        t = surv_df[, W], 
                                        surv = surv_df[, "surv0"], 
                                        surv_between = surv_between)
  if (surv_beyond == "w") {
    SURVmax = data[max(which(uncens)), "surv0"]
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
      data$imp = NA
      return(list(imputed_data = data, code = FALSE))
    }
    needs_extrap = which(!uncens & data[, W] > Xtilde)
    data[needs_extrap, "surv0"] = sapply(X = data[needs_extrap, W], 
                                         FUN = extrap_surv_beyond, 
                                         t = surv_df[, W], 
                                         surv = surv_df[, "surv0"], 
                                         surv_beyond = surv_beyond, 
                                         weibull_params = weibull_params)
  } else {
    needs_extrap = which(!uncens & data[, W] > Xtilde)
    data[needs_extrap, "surv0"] = sapply(X = data[needs_extrap, W], 
                                         FUN = extrap_surv_beyond, 
                                         t = surv_df[, W], 
                                         surv = surv_df[, "surv0"], 
                                         surv_beyond = surv_beyond)
    weibull_params = NULL
  }
  data$surv = data[, "surv0"] ^ data[, "HR"]
  if (trapezoidal_rule) {
    data_dist = unique(data[, c(W, Delta, Z, "surv0")])
    t_diff = data_dist[-1, W] - data_dist[-nrow(data_dist), W]
    for (i in which(!uncens)) {
      surv_sum_i = data_dist[-1, "surv0"] ^ as.numeric(data[i, "HR"]) + 
        data_dist[-nrow(data_dist), "surv0"] ^ as.numeric(data[i, "HR"]) 
      int_surv_i = 1 / 2 * sum((data_dist[-nrow(data_dist), W] >= as.numeric(data[i, W])) * surv_sum_i * t_diff)
      data$imp[i] = data$imp[i] + (int_surv_i / data[i, "surv"])
    }
  } else {
    to_integrate = function(t, hr) {
      basesurv = sapply(X = t, 
                        FUN = extend_surv, 
                        t = surv_df[, W], 
                        surv = surv_df[, "surv0"], 
                        surv_between = surv_between, 
                        surv_beyond = surv_beyond, 
                        weibull_params = weibull_params)
      basesurv ^ as.numeric(hr)
    }
    int_surv = sapply(X = which(!uncens), 
                      FUN = function(i) {
                        tryCatch(expr = integrate(f = to_integrate, 
                                                  lower = data[i, W], 
                                                  upper = Xmax, 
                                                  subdivisions = subdivisions, 
                                                  hr = data[i, "HR"], 
                                                  rel.tol = .Machine$double.eps ^ 0.24)$value, 
                                 error = function(e) return(NA))}
    )
    data$imp[which(!uncens)] = data[which(!uncens), W] + int_surv / data[which(!uncens), "surv"]
  }
  
  return(list(imputed_data = data, 
              code = !any(is.na(data$imp))))
}