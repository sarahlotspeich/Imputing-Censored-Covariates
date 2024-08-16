#' Semiparametric conditional mean imputation (CMI) for a censored predictor
#'
#' Semiparametric conditional mean imputation (CMI) for a censored predictor using a Cox proportional hazards model and the Breslow estimator to estimate conditional survival.
#'
#' @param imputation_model Imputation model formula (or coercible to formula), a formula expression as for other regression models. The response is usually a survival object as returned by the \code{Surv} function. See the documentation for \code{Surv} for details.
#' @param lp (Optional) A vector of the linear predictors for the imputation model, in the same order as \code{data}. If \code{lp = NULL} (the default), the linear predictor is extracted from the fitted imputation model.
#' @param data Dataframe or named matrix containing columns \code{W}, \code{Delta}, and \code{Z}.
#' @param integral A string input for how to approximate the integral in the imputed values. Default is \code{integral="aq"} for adaptive quadrature, but \code{"tr"} (trapezoidal rule) and \code{"a"} (quasi-analytical) are also available.
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
#' @importFrom expint gammainc

cmi_sp = function (imputation_model, lp = NULL, data, integral = "aq", Xmax = Inf, subdivisions = 2000, surv_between = "cf", surv_beyond = "e") {
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
  
  if (is.null(lp)) {
    # Calculate linear predictor for Cox imputation model
    lp = predict(fit, reference = "sample") + sum(coef(fit) * fit$means, na.rm = TRUE) ## linear predictors
  }
  
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
  data$imp = data[, W]
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
  if (integral == "tr") {
    data_dist = unique(data[, c(W, Delta, Z, "surv0")])
    t_diff = data_dist[-1, W] - data_dist[-nrow(data_dist), W]
    for (i in which(!uncens)) {
      surv_sum_i = data_dist[-1, "surv0"] ^ as.numeric(data[i, "HR"]) + 
        data_dist[-nrow(data_dist), "surv0"] ^ as.numeric(data[i, "HR"]) 
      int_surv_i = 1 / 2 * sum((data_dist[-nrow(data_dist), W] >= as.numeric(data[i, W])) * surv_sum_i * t_diff)
      data$imp[i] = data$imp[i] + (int_surv_i / data[i, "surv"])
    }
  } else if (integral == "aq") {
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
                                                   hr = data[i, "HR"])$value, 
                                  error = function(e) return(NA))}
    )
    data$imp[which(!uncens)] = data[which(!uncens), W] + int_surv/data[which(!uncens), "surv"]
  } else if (integral == "a") {
    # Estimate the integral up to Xtilde using the trapezoidal rule 
    data_dist = unique(data[uncens, c(W, Delta, Z, "surv0")])
    t_diff = data_dist[-1, W] - data_dist[-nrow(data_dist), W]
    tr = vector()
    for (i in which(!uncens)) {
      surv_sum_i = data_dist[-1, "surv0"] ^ as.numeric(data[i, "HR"]) +
        data_dist[-nrow(data_dist), "surv0"] ^ as.numeric(data[i, "HR"])
      tr = append(tr,
                  1 / 2 * sum((data_dist[-nrow(data_dist), W] >= as.numeric(data[i, W])) * surv_sum_i * t_diff))
    }
    # to_integrate = function(t, hr) {
    #   basesurv = sapply(X = t, 
    #                     FUN = extend_surv, 
    #                     t = surv_df[, W], 
    #                     surv = surv_df[, "surv0"], 
    #                     surv_between = surv_between, 
    #                     surv_beyond = surv_beyond, 
    #                     weibull_params = weibull_params)
    #   basesurv ^ as.numeric(hr)
    # }
    # aq = sapply(X = which(!uncens), 
    #             FUN = function(i) {
    #               tryCatch(expr = integrate(f = to_integrate, 
    #                                         lower = data[i, W], 
    #                                         upper = Xtilde, 
    #                                         subdivisions = subdivisions, 
    #                                         hr = data[i, "HR"])$value, 
    #                        error = function(e) return(NA))}
    #             )
    
    # Estimate the integral from Xtilde to infinity using the trapezoidal rule 
    if (surv_beyond == "w") {
      ## Get parameter estimates 
      alpha_hat = weibull_params[1]
      lambda_tilde = weibull_params[2] * data[which(!uncens), "HR"]
      
      ## Get integral estimates 
      ### Integral from Xtilde to infinity
      a1 = gammainc(a = 1 / alpha_hat, x = lambda_tilde * Xtilde ^ alpha_hat) / 
        (lambda_tilde * Xtilde ^ (1 / alpha_hat))
      ### Integral from Wi to infinity
      a2 = gammainc(a = 1 / alpha_hat, x = lambda_tilde * data[which(!uncens), "w"] ^ alpha_hat) / 
        (lambda_tilde * data[which(!uncens), "w"] ^ (1 / alpha_hat))
      ### If Wi <= Xtilde, take a1; If Wi > Xtilde; take a2
      a = pmin(a1, a2)
    } else if (surv_beyond == "e") {
      ## Get parameter estimates 
      rho_hat = - Xtilde / log(SURVmax)
      rho_tilde = data[which(!uncens), "HR"] / rho_hat 
      
      ## Get integral estimates 
      ### Integral from Xtilde to infinity
      a1 = exp(- rho_tilde * Xtilde) / rho_tilde
      ### Integral from Wi to infinity
      a2 = exp(- rho_tilde * data[which(!uncens), "w"]) / rho_tilde
      ### If Wi <= Xtilde, take a1; If Wi > Xtilde; take a2
      a = pmin(a1, a2)
    }
    
    # Take sum of integrals
    # int_surv = aq + a 
    int_surv = tr + a
    
    # Compute conditional means
    data$imp[which(!uncens)] = data[which(!uncens), W] + int_surv / data[which(!uncens), "surv"]
  }
  if (any(is.na(data$imp))) {
    data$imp[which(is.na(data$imp))] = data[which(is.na(data$imp)), W]
  }
  if (any(data$imp == Inf)) {
    data$imp[which(data$imp == Inf)] = data[which(data$imp == Inf), W]
  }
  
  return(list(imputed_data = data, code = !any(is.na(data$imp))))
}