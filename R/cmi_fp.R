#' Fully parametric conditional mean imputation (CMI) for a censored predictor
#'
#' Fully parametric conditional mean imputation (CMI) for a censored predictor using an accelerated failure-time (AFT) model to estimate conditional survival.
#'
#' @param W character, column name for observed values of the censored covariate 
#' @param Delta character, column name for censoring indicators. Note that \code{Delta = 0} is interpreted as a censored observation. 
#' @param Z character vector, column name(s) of additional fully observed covariate(s).
#' @param data Dataframe or named matrix containing columns \code{W}, \code{Delta}, and \code{Z}.
#' @param fit \code{survreg} fit, imputation model object modeling \code{W} on \code{Z}. If \code{fit = NULL} (default), the AFT model with only main effects for \code{Z} and assuming distribution \code{dist} is fit internally and used.
#' @param dist (optional) character, assumed distribution for \code{W} in the AFT model, passed to \code{survival::survreg()}. If \code{fit} is supplied, no need to supply \code{dist}; if \code{fit} is not supplied, default is \code{dist = "weibull"}.
#' @param infinite_integral (optional) logical, if \code{infinite_integral = TRUE} (default) then conditional means are found by integrating from \code{W} to \code{Inf}
#' @param max_iter (optional) numeric, maximum iterations allowed in call to \code{survival::survreg()}. Default is \code{100}.
#'
#' @return 
#' \item{imputed_data}{A copy of \code{data} with added column \code{imp} containing the imputed values.}
#' \item{code}{Indicator of algorithm status (\code{TRUE} or \code{FALSE}).}
#'
#' @export
#' @importFrom survival survreg 
#' @importFrom survival Surv 
#' @importFrom survival psurvreg

cmi_fp = function(W, Delta, Z, data, fit = NULL, dist = "weibull", 
                  infinite_integral = TRUE, maxiter = 100) {
  # If no imputation model was supplied, fit an AFT model using main effects only
  if (is.null(fit)) {
    fit_formula = as.formula(paste0("Surv(time = ", W, ", event = ", Delta, ") ~ ", paste0(Z, collapse = " + ")))
    fit = tryCatch(expr = survreg(formula = fit_formula, 
                                  data = data, 
                                  dist = dist, 
                                  control = list(maxiter = maxiter)),
                   warning = function(w) return(NA))
  }
  
  # If the imputation model does not converge, we cannot impute 
  if (any(is.na(fit))) {
    return(list(imputed_data = data, code = FALSE))
  }
  
  # Calculate linear predictor \lambda %*% Z for AFT model
  lp = fit$coefficients[1] +
    data.matrix(data[, Z]) %*% matrix(data = fit$coefficients[- 1], ncol = 1)
 
  # Calculate survival with original model coefficients using built-in function psurvreg (returns CDF)
  data = data.frame(data, 
                    surv = 1 - psurvreg(q = data[, W], 
                                        mean = lp, 
                                        scale = fit$scale, 
                                        distribution = dist))

  # Order data by W
  data = data[order(data[, W]), ]
  
  # Create an indicator variable for being uncensored
  uncens = data[, Delta] == 1
  
  # Calculate imputed values 
  data$imp = data[, W]
  
  if (infinite_integral) {
    ## Use integrate() to approximate integral from W to \infty of S(t|Z)
    int_surv = sapply(
      X = which(!uncens), 
      FUN = function(i) { 
        tryCatch(expr = integrate(f = function(t) { 1 - psurvreg(q = t, mean = lp[i], scale = fit$scale, distribution = dist) }, 
                                  lower = data[i, W], 
                                  upper = Inf)$value,
                 error = function(e) return(NA))
      }
    )
  } else {
    ## Calculate mean life = integral from 0 to \infty of S(t|Z)
    if (dist %in% c("weibull", "exponential", "rayleigh")) {
      est_ml = exp(lp[!uncens]) * gamma(1 + fit$scale)
    } #else if (dist %in% c("lognormal", "loggaussian")) {}
    
    ## Use integrate() to approximate integral from 0 to W of S(t|Z)
    int_surv = sapply(
      X = which(!uncens), 
      FUN = function(i) { 
        tryCatch(expr = integrate(function(t) { 1 - psurvreg(q = t, mean = lp[i], scale = fit$scale, distribution = dist) }, 
                                  lower = 0, 
                                  upper = data[i, W])$value,
                 error = function(e) return(NA))
      }
    )
    
    ## Subtract integral from mean life to get integral from W to \infty of S(t|Z)
    int_surv = est_ml - int_surv
  }
  
  ## Calculate MRL = int_surv / surv(W|Z)
  data$mrl = 0
  data$mrl[which(!uncens)] = int_surv / data[which(!uncens), "surv"]
  
  ## Calculate E(X|X>W,Z) = W + int_surv / surv(W|Z)
  data$imp[which(!uncens)] = data[which(!uncens), W] + data$mrl[which(!uncens)]
  
  ## Check for infinite/NA imputed values 
  if (any(is.na(data$imp))) {
    data$imp[which(is.na(data$imp))] = data[which(is.na(data$imp)), W]
  }
  if (any(data$imp == Inf)) {
    data$imp[which(data$imp == Inf)] = data[which(data$imp == Inf), W]
  }
  
  # Return input dataset with appended column imp containing imputed values 
  if (any(is.na(data$imp))) {
    return(list(imputed_data = data, code = FALSE))  
  } else {
    return(list(imputed_data = data, code = TRUE))
  }
}