#' Semiparametric conditional mean imputation (CMI) for a censored predictor
#'
#' Semiparametric conditional mean imputation (CMI) for a censored predictor using a Cox proportional hazards model and the Breslow estimator to estimate conditional survival. 
#' Regression coefficients from the Cox model are resampled from their estimated multivariate normal distribution.
#'
#' @param imputation_model Imputation model formula (or coercible to formula), a formula expression as for other regression models. The response is usually a survival object as returned by the \code{Surv} function. See the documentation for \code{Surv} for details.#' @param data Dataframe or named matrix containing columns \code{W}, \code{Delta}, and \code{Z}.
#' @param M Number of times to resample Cox model coefficients
#' @param trapezoidal_rule A logical input for whether the trapezoidal rule should be used to approximate the integral in the imputed values. Default is \code{FALSE}.
#' @param Xmax (Optional) Upper limit of the domain of the censored predictor. Default is \code{Xmax = Inf}.
#' @param surv_between A string for the method to be used to interpolate for censored values between events. Options include \code{"cf"} (carry forward, the default), \code{"wm"} (weighted mean), or \code{"m"} (mean).
#' @param surv_beyond A string for the method to be used to extrapolate the survival curve beyond the last observed event. Options include \code{"d"} (immediate drop off), \code{"e"} (exponential extension, the default), or \code{"w"} (weibull extension).
#' 
#' @return 
#' \item{imputed_datasets}{An list of length m, containing copies of \code{data} each with added column \code{imp} containing the imputed values.}
#'
#' @export
#' @importFrom survival coxph 
#' @importFrom survival Surv
#' @importFrom MASS mvrnorm
cmi_sp_resample = function (imputation_model, data, M = 1, trapezoidal_rule = FALSE, Xmax = Inf, surv_between = "cf", surv_beyond = "e") {
  # Extract variable names from imputation_model
  W = all.vars(imputation_model)[1] ## censored covariate
  Delta = all.vars(imputation_model)[2] ## corresponding event indicator
  Z = all.vars(imputation_model)[-c(1:2)] ## additional covariates
  
  # Fit Cox PH imputation model for X ~ Z 
  fit = coxph(formula = imputation_model, 
              data = data)
  
  # resample Cox coefficients from estimated normal distribution
  if (M == 1) {
    Cox_coef = as.matrix(coef(fit), ncol = 1)
  } else {
    Cox_coef = t(mvrnorm(n = M, 
                         mu = coef(fit), 
                         Sigma = vcov(fit)))
    }
  
  # store imputed datasets in list
  imputed_datasets = list()
  
  for (m in 1:M) {
    lp = as.matrix(x = data[, Z]) %*% Cox_coef[, m]
    
    data$HR = exp(lp)
    be = breslow_estimator(x = NULL, 
                           time = W, 
                           event = Delta, 
                           hr = "HR", 
                           data = data)
    surv_df = with(be, 
                   data.frame(t = times, surv0 = basesurv))
    colnames(surv_df)[1] <- W
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
      return(list(imputed_data = data, 
                  code = FALSE))
    }
    needs_extrap = which(!uncens & data[, W] > Xtilde)
    data[needs_extrap, "surv0"] <- sapply(X = data[needs_extrap, W], 
                                          FUN = extrap_surv_beyond, 
                                          t = surv_df[, W], 
                                          surv = surv_df[, "surv0"], 
                                          surv_beyond = surv_beyond, 
                                          weibull_params = weibull_params)
    } else {
      needs_extrap <- which(!uncens & data[, W] > Xtilde)
      data[needs_extrap, "surv0"] <- sapply(X = data[needs_extrap, W], 
                                            FUN = extrap_surv_beyond, 
                                            t = surv_df[, W], 
                                            surv = surv_df[, "surv0"], 
                                            surv_beyond = surv_beyond)
      weibull_params <- NULL
    }
    data$surv = data[, "surv0"] ^ data[, "HR"]
    if (trapezoidal_rule) {
      data_dist = unique(data[, c(W, Delta, Z, "surv")])
      t_diff = data_dist[-1, W] - data_dist[-nrow(data_dist), W]
      surv_sum = data_dist[-1, "surv"] + data_dist[-nrow(data_dist), "surv"]
      for (i in which(!uncens)) {
        sum_surv_i = sum((data_dist[-nrow(data_dist), W] >= as.numeric(data[i, W])) * surv_sum * t_diff)
        data$imp[i] = data$imp[i] + (1/2) * (sum_surv_i/data[i, "surv"])
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
                                                    subdivisions = 2000, 
                                                    hr = data[i, "HR"])$value, 
                                    error = function(e) return(NA))}
      )
      data$imp[which(!uncens)] = data[which(!uncens), W] + int_surv/data[which(!uncens), "surv"]
    }
    if (any(is.na(data$imp))) {
      data$imp[which(is.na(data$imp))] = data[which(is.na(data$imp)), W]
    }
    if (any(data$imp == Inf)) {
      data$imp[which(data$imp == Inf)] = data[which(data$imp == Inf), W]
    }
    
    imputed_datasets[[m]] = data
  }
  
  return(imputed_datasets)
}