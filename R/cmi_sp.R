#' Semiparametric conditional mean imputation (CMI) for a censored predictor
#'
#' Semiparametric conditional mean imputation (CMI) for a censored predictor using a Cox proportional hazards model and the Breslow estimator to estimate conditional survival.
#'
#' @param W Column name of observed predictor values (including censored opens). 
#' @param Delta Column name of censoring indicators. Note that \code{data[, Delta] = 0} is interpreted as a censored observation. 
#' @param Z Column name of additional fully observed covariates.
#' @param data Dataframe or named matrix containing columns \code{W}, \code{Delta}, and \code{Z}.
#' @param fit (Optional) A \code{coxph} imputation model object modeling \code{W} on \code{Z}. If \code{fit = NULL} (default), the Cox model with only main effects for \code{Z} is used.
#' @param stratified If \code{TRUE}, stratification in \code{W} is used to construct time-varying coefficients. Default is \code{FALSE}. 
#' @param split_data (Optional) If \code{stratified = TRUE} and \code{fit} is supplied, additional dataframe of stratified data.
#' @param trapezoidal_rule A logical input for whether the trapezoidal rule should be used to approximate the integral in the imputed values. Default is \code{FALSE}.
#' @param surv_between A string for the method to be used to interpolate for censored values between events. Options include \code{"carry-forward"} (default), \code{"linear"}, or \code{"mean"}.
#' @param surv_beyond A string for the method to be used to extrapolate the survival curve beyond the last observed event. Options include \code{"drop-off"}, \code{"exponential"} (default), or \code{"weibull"}.
#'
#' @return 
#' \item{imputed_data}{A copy of \code{data} with added column \code{imp} containing the imputed values.}
#' \item{code}{Indicator of algorithm status (\code{TRUE} or \code{FALSE}).}
#' \item{splits}{(If \code{stratified = TRUE}) The number of strata used.}
#'
#' @export
#' @importFrom survival coxph 
#' @importFrom survival Surv
#' @importFrom survival survSplit
#' @importFrom survival cox.zph

cmi_sp <- function(W, Delta, Z, data, fit = NULL, stratified = FALSE, split_data = NULL, trapezoidal_rule = FALSE, surv_between = "carry-forward", surv_beyond = "exponential") {
  fit_internally <- is.null(fit)
  
  # If no imputation model was supplied, fit a Cox PH using main effects
  if (fit_internally) {
    fit_formula <- as.formula(paste0("Surv(time = ", W, ", event = ", Delta, ") ~ ", paste0(Z, collapse = " + ")))
    fit <- coxph(formula = fit_formula, 
                 data = data)
  }
  
  if (stratified & fit_internally) {
    # Find the smallest number of splits needed for PH to be upheld
    splits <- 1
    testPH <- capture.output(print(cox.zph(fit = fit))[1, 3])
    pval <- as.numeric(substr(x = testPH[4], start = 5, stop = nchar(testPH[4])))
    while (pval < 0.05) {
      # Split data into percentile-based intervals
      splits <- splits + 1
      where_split <- seq(from = 0, to = 1, by = 1 / splits)
      where_split <- where_split[-c(1, length(where_split))]
      split_dat <- survSplit(formula = fit_formula, data = data,
                             cut = quantile(data[, W], where_split),
                             episode = "tgroup", id = "id")
      # Fit the Cox model to the split data
      split_fit_formula <- as.formula(paste0("Surv(time = tstart, time2 = ", W, ", event = ", Delta, ") ~ ", paste0(Z, collapse = " + "), ":strata(tgroup)"))
      fit <- coxph(formula = split_fit_formula, data = split_dat, control = coxph.control(timefix = FALSE))
      testPH <- capture.output(print(cox.zph(fit = fit))[1, 3])
      pval <- as.numeric(substr(x = testPH[4], start = 5, stop = nchar(testPH[4])))
    }
  } 
  
  # Calculate linear predictor \lambda %*% Z for Cox model
  if (stratified) {
    #if (splits > 1) {}
    suppressWarnings(
      split_lp <- predict(fit, reference="sample") + sum(coef(fit) * fit$means, na.rm = TRUE)
    )
    agg_lp <- aggregate(x = split_lp ~ split_dat$id, FUN = sum)
    lp <- agg_lp[, 2]
  } else {
    lp <- predict(fit, reference="sample") + sum(coef(fit) * fit$means, na.rm = TRUE)
  }
  
  ## Calculate hazard ratio for Cox model
  data$HR <- exp(lp)
  
  # Estimate baseline survival from Cox model fit using Breslow's estimator
  be <- breslow_estimator(x = NULL, 
                          time = W, 
                          event = Delta, 
                          hr = "HR", 
                          data = data)
  surv_df <- with(be, data.frame(t = times, surv0 = basesurv))
  colnames(surv_df)[1] <- W
  
  ## Merge baseline survival estimates into data
  data <- merge(x = data, y = surv_df, all.x = TRUE, sort = FALSE)
  
  # Order data by W
  data <- data[order(data[, W]), ]
  
  # Create an indicator variable for being uncensored
  uncens <- data[, Delta] == 1
  
  # For people with events, obs = X
  data$imp <- data[, W]
  
  # Assume survival at censored W < min(X) = 1
  minX <- data[min(which(uncens)), W] 
  data[which(data[, W] < minX), "surv0"] <- 1
  
  # Interpolate baseline survival at censored W < \widetilde{X}
  Xtilde <- data[max(which(uncens)), W] 
  needs_interp <- which(is.na(data[, "surv0"]) & data[, W] < Xtilde)
  data[needs_interp, "surv0"] <- sapply(X = data[needs_interp, W], 
                                        FUN = interp_surv_between, 
                                        t = surv_df[, W], 
                                        surv = surv_df[, "surv0"], 
                                        surv_between = surv_between)
  
  # Extrapolate baseline survival at censored W > \widetilde{X}
  if (surv_beyond == "weibull") {
    # Estimate Weibull parameters using constrained MLE
    SURVmax <- data[max(which(uncens)), "surv0"]
    weibull_params <- constr_weibull_mle(t = data[, W], I_event = data[, Delta], Xtilde = Xtilde, rho = SURVmax, alpha0 = 1E-4)
    
    # If weibull params don't converge, quit 
    if (any(is.na(weibull_params))) {
      data$imp <- NA
      return(list(imputed_data = data, code = FALSE))   
    }
    
    # If they do, extrapolate with them 
    needs_extrap <- which(!uncens & data[, W] > Xtilde)
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

  # Calculate conditional survival W | Z for all 
  data$surv <- data[, "surv0"] ^ data[, "HR"]
  
  # Calculate imputed values E(X|X>W,Z)
  if (trapezoidal_rule) {
    # Distinct rows (in case of non-unique obs values)
    data_dist <- unique(data[, c(W, Delta, Z, "surv")])
    
    # [T_{(i+1)} - T_{(i)}]
    t_diff <- data_dist[- 1, W] - data_dist[- nrow(data_dist), W]
    
    # S(t+1|Z) + S(t|Z)
    surv_sum <- data_dist[-1, "surv"] + data_dist[- nrow(data_dist), "surv"]
    
    # Use trapezoidal approximation for integral
    for (i in which(!uncens)) {
      sum_surv_i <- sum((data_dist[-nrow(data_dist), W] >= as.numeric(data[i, W])) * surv_sum * t_diff)
      data$imp[i] <- data$imp[i] + (1 / 2) * (sum_surv_i /  data[i, "surv"])
    }
  } else {
    # Builds on the extend_surv function by raising S_0(t) ^ HR = S(t|Z)
    to_integrate <- function(t, hr) {
      basesurv <- sapply(X = t, 
                         FUN = extend_surv, 
                         t = surv_df[, W], 
                         surv = surv_df[, "surv0"], 
                         surv_between = surv_between, 
                         surv_beyond = surv_beyond, 
                         weibull_params = weibull_params)  
      basesurv ^ as.numeric(hr)
    }
    
    ## Use integrate() to approximate integral from W to \infty of S(t|Z)
    int_surv <- sapply(
      X = which(!uncens), 
      FUN = function(i) { 
        tryCatch(expr = integrate(f = to_integrate, 
                                  lower = data[i, W], 
                                  upper = Inf, 
                                  subdivisions = 2000, 
                                  hr = data[i, "HR"])$value,
                 error = function(e) return(NA))
      }
    )
    
    ## Calculate E(X|X>W,Z) = W + int_surv / surv(W|Z)
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
    if (stratified) {
      return(list(imputed_data = data, code = FALSE, splits = splits))  
    } else {
      return(list(imputed_data = data, code = FALSE, splits = NA))    
    }
  } else {
    if (stratified) {
      return(list(imputed_data = data, code = TRUE, splits = splits))  
    } else {
      return(list(imputed_data = data, code = TRUE, splits = NA))
    }
  }
}