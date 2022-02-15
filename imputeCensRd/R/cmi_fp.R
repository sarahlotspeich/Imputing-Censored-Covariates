#' Fully parametric conditional mean imputation (CMI) for a censored predictor
#'
#' Fully parametric conditional mean imputation (CMI) for a censored predictor using an accelerated failure-time (AFT) model to estimate conditional survival.
#'
#' @param W Column name of observed predictor values (including censored opens). 
#' @param Delta Column name of censoring indicators. Note that \code{Delta = 0} is interpreted as a censored observation. 
#' @param Z Column name of additional fully observed covariates.
#' @param data Dataframe or named matrix containing columns \code{W}, \code{Delta}, and \code{Z}.
#' @param fit A \code{survreg} imputation model object modeling \code{W} on \code{Z}. If \code{fit = NULL} (default), the AFT model with only main effects for \code{Z} and assuming a Weibull distribution is fit internally and used.
#' @param dist (Optional) Assumed distribution for \code{W} in the AFT model, passed to \code{survival::survreg()}. Default is \code{"weibull"}.
#'
#' @return 
#' \item{imputed_data}{A copy of \code{data} with added column \code{imp} containing the imputed values.}
#' \item{code}{Indicator of algorithm status (\code{TRUE} or \code{FALSE}).}
#'
#' @export
#' @importFrom survival survreg Surv psurvreg

cmi_fp <- function(W, Delta, Z, data, fit = NULL, dist = "weibull") {
  # If no imputation model was supplied, fit an AFT model using main effects
  if (is.null(fit)) {
    fit_formula <- as.formula(paste0("Surv(time = ", W, ", event = ", Delta, ") ~ ", paste0(Z, collapse = " + ")))
    fit <- survreg(formula = Surv(time = t, event = d) ~ z, data = data, dist = dist)
  }
  
  # If the imputation model does not converge, we cannot impute 
  if (any(is.na(fit$coef))) {
    return(list(imputed_data = data, code = FALSE))
  }
  
  # Calculate linear predictor \lambda %*% Z for AFT model
  lp <- fit$coefficients[1] + 
    data.matrix(data[, Z]) %*% matrix(data = fit$coefficients[-1], ncol = 1)
 
  # Calculate survival with original model coefficients using built-in function
  surv_df <- data.frame(t = data[, W], surv = 1 - psurvreg(q = data[, W], mean = lp, scale = fit$scale, distribution = dist))
  colnames(surv_df)[1] <- W
  
  # Merge survival estimates into data
  data <- merge(x = data, y = surv_df, all.x = TRUE, sort = FALSE)
  
  # Order data by W
  data <- data[order(data[, W]), ]
  
  # Create an indicator variable for being uncensored
  uncens <- data[, Delta] == 1
  
  # Calculate E(X|X>C,Z) using integrate() 
  data$imp <- data[, W]
  data$imp[which(!uncens)] <- sapply(X = which(!uncens), 
                               FUN = function(i) { 
                                 tryCatch(expr = integrate(f = function(t) 1 - psurvreg(q = t, mean = lp[i], scale = fit$scale, distribution = dist), lower = data[i, W], upper = Inf)$value,
                                          error = function(e) return(NA))
                               })
  
  # Return input dataset with appended column imp containing imputed values 
  if (any(is.na(data$imp))) {
    return(list(imputed_data = data, code = FALSE))  
  } else {
    return(list(imputed_data = data, code = TRUE))
  }
}