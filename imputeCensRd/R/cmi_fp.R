#' Fully parametric conditional mean imputation (CMI) for a censored predictor
#'
#' Fully parametric conditional mean imputation (CMI) for a censored predictor using an accelerated failure-time (AFT) model to estimate conditional survival.
#'
#' @param W Column name of observed predictor values (including censored opens). 
#' @param Delta Column name of censoring indicators. Note that \code{Delta = 0} is interpreted as a censored observation. 
#' @param Z Column name of additional fully observed covariates.
#' @param data Dataframe or named matrix containing columns \code{W}, \code{Delta}, and \code{Z}.
#' @param fit A \code{survreg} imputation model object modeling \code{W} on \code{Z}. If \code{fit = NULL} (default), the AFT model with only main effects for \code{Z} and assuming a Weibull distribution is fit internally and used.
#' @param dist Assumed distribution for \code{W} in the AFT model, passed to \code{survival::survreg()}. Default is \code{"weibull"}.
#' @param trapezoidal_rule A logical input for whether the trapezoidal rule should be used to approximate the integral in the imputed values. Default is \code{FALSE}.
#' @param max_iter Maximum iterations allowed in call to \code{survival::survreg()}. Default is \code{100}.
#'
#' @return 
#' \item{imputed_data}{A copy of \code{data} with added column \code{imp} containing the imputed values.}
#' \item{code}{Indicator of algorithm status (\code{TRUE} or \code{FALSE}).}
#'
#' @export
#' @importFrom survival survreg 
#' @importFrom survival Surv 
#' @importFrom survival psurvreg

cmi_fp <- function(W, Delta, Z, data, fit = NULL, dist = "weibull", trapezoidal_rule = FALSE, maxiter = 100) {
  # If no imputation model was supplied, fit an AFT model using main effects
  if (is.null(fit)) {
    fit_formula <- as.formula(paste0("Surv(time = ", W, ", event = ", Delta, ") ~ ", paste0(Z, collapse = " + ")))
    fit <- tryCatch(expr = survreg(formula = fit_formula, data = data, dist = dist, control = list(maxiter = maxiter)),
                    warning = function(w) return(NA))
  }
  
  # If the imputation model does not converge, we cannot impute 
  if (any(is.na(fit))) {
    return(list(imputed_data = data, code = FALSE))
  }
  
  # Calculate linear predictor \lambda %*% Z for AFT model
  lp <- fit$coefficients[1] + 
    data.matrix(data[, Z]) %*% matrix(data = fit$coefficients[- 1], ncol = 1)
 
  # Calculate survival with original model coefficients using built-in function
  data <- data.frame(data, surv = 1 - psurvreg(q = data[, W], mean = lp, scale = fit$scale, distribution = dist))

  # Order data by W
  data <- data[order(data[, W]), ]
  
  # Create an indicator variable for being uncensored
  uncens <- data[, Delta] == 1
  
  # Calculate imputed values 
  data$imp <- data[, W]
  if (trapezoidal_rule) {
    # Distinct rows (in case of non-unique obs values)
    data_dist <- unique(data[, c(W, Delta, Z, "surv")])
    
    # [T_{(i+1)} - T_{(i)}]
    t_diff <- data_dist[- 1, W] - data_dist[- nrow(data_dist), W]
    
    # Censored subject values (to impute)
    t_cens <- data[data[, Delta] == 0, W]
    
    # Follow formula assuming AFT model for S(X|Z)
    for (x in which(!uncens)) {
      Cj <- data[x, W]
      Sj <- data_dist[-1, "surv"] + data_dist[- nrow(data_dist), "surv"]
      num <- sum((data_dist[-nrow(data_dist), W] >= Cj) * Sj * t_diff)
      denom <- data[x, "surv"]
      data$imp[x] <- (1 / 2) * (num / denom) + Cj
    }
  } else {
    ## Use integrate() to approximate integral from W to \infty of S(t|Z)
    int_surv <- sapply(
      X = which(!uncens), 
      FUN = function(i) { 
        tryCatch(expr = integrate(f = function(t) 1 - psurvreg(q = t, mean = lp[i], scale = fit$scale, distribution = dist), lower = data[i, W], upper = Inf)$value,
                 error = function(e) return(NA))
      }
    )
    ## Calculate E(X|X>W,Z) = int_surv / surv(W|Z) + W
    data$imp[which(!uncens)] <- data[which(!uncens), W] + int_surv / data[which(!uncens), "surv"]
  }
  
  ## Check for infinite/NA imputed values 
  if (any(is.na(data$imp))) {
    data$imp[which(is.na(data$imp))] <- data[which(is.na(data$imp)), W]
  }
  if (any(data$imp  == Inf)) {
    data$imp[which(data$imp == Inf)] <- data[which(data$imp == Inf), W]
  }
  
  # Return input dataset with appended column imp containing imputed values 
  if (any(is.na(data$imp))) {
    return(list(imputed_data = data, code = FALSE))  
  } else {
    return(list(imputed_data = data, code = TRUE))
  }
}