#' Semiparametric conditional mean imputation (CMI) for a censored predictor
#'
#' Semiparametric conditional mean imputation (CMI) for a censored predictor using a Cox proportional hazards model and the Breslow estimator to estimate conditional survival.
#'
#' @param W Column name of observed predictor values (including censored opens). 
#' @param Delta Column name of censoring indicators. Note that \code{data[, Delta] = 0} is interpreted as a censored observation. 
#' @param Z Column name of additional fully observed covariates.
#' @param data Dataframe or named matrix containing columns \code{W}, \code{Delta}, and \code{Z}.
#' @param fit A \code{coxph} imputation model object modeling \code{W} on \code{Z}. If \code{fit = NULL} (default), the Cox model with only main effects for \code{Z} is fit internally and used.
#' @param extrapolate A string for the method to be used to extrpolate the survival curve beyond the last observed event. Options include \code{"none"} (default), \code{"exponential"}, or \code{"weibull"}.
#' @param interpolate_between A string for the method to be used to interpolate for censored values between events. Options include \code{"carry-forward"} (default), \code{"linear"}, or \code{"mean"}.
#' @param forceLastEvent A logical input to force the last observed value for \code{W} to be treated as an event. Default is \code{FALSE}.
#'
#' @return A copy of \code{data} with added column \code{imp} containing the imputed values.
#'
#' @export
#' @importFrom survival coxph Surv survreg psurvreg

cmi_sp <- function(W, Delta, Z, data, fit = NULL, extrapolate = "none", interpolate_between = "carry-forward", forceLastEvent = FALSE) {
  # Assume last observed value is an event regardless
  if (forceLastEvent) {
    data[which.max(data[, W]), Delta] <- 1
  }
  
  # If no imputation model was supplied, fit a Cox PH using main effects
  if (is.null(fit)) {
    fit_formula <- as.formula(paste0("Surv(time = ", W, ", event = ", Delta, ") ~ ", paste0(Z, collapse = " + ")))
    fit <- coxph(formula = fit_formula, data = data)
  }
  
  # Calculate linear predictor \lambda %*% Z for Cox model
  lp <- data.matrix(data[, Z]) %*% matrix(data = fit$coefficients, ncol = 1)
  data$HR <- exp(lp)
  
  # Estimate baseline survival from Cox model fit using Breslow estimator
  be <- breslow_estimator(x = NULL, time = W, event = Delta, hr = "HR", data = data)
  surv_df <- with(be, data.frame(t = times, surv0 = basesurv))
  colnames(surv_df)[1] <- W
  
  # Merge baseline survival estimates into data
  data <- merge(x = data, y = surv_df, all.x = TRUE, sort = FALSE)
  
  # Calculate conditional survival
  data$surv <- with(data, surv0 ^ HR)
  
  # Order data by W
  data <- data[order(data[, W]), ]
  
  # Create an indicator variable for being uncensored
  uncens <- data[, Delta] == 1
  
  # For people with events, obs = X
  data$imp <- data[, W]
  
  if (extrapolate == "none") {
    if (any(is.na(data[, "surv"]))) {
      # For censored subjects, survival is average of times right before/after
      suppressWarnings(
        data[is.na(data[, "surv"]), "surv"] <- sapply(X = data[is.na(data[, "surv"]), W], FUN = impute_censored_surv, time = W, event = Delta, surv = "surv", data = data)
      )
    }
    
    # Assume survival = 1 before first observed covariate value
    if (any(data[, W] < min(data[uncens, W]))) {
      cens_before <- which(data[, W] < min(data[uncens, W]))
      data[cens_before, "surv"] <- 1
    }
    
    # Extrapolate survival beyond last observed covariate
    if (any(data[, W] > max(data[uncens, W]))) {
      cens_after <- which(data[, W] > max(data[uncens, W]))
      t_cens_after <- data[cens_after, W]
      last_event_surv <- data[max(which(uncens)), "surv"]
      # Gill (1980) carry forward survival at last event
      data[cens_after, "surv"] <- last_event_surv
    }
    
    # Distinct rows (in case of non-unique obs values)
    data_dist <- unique(data[, c(W, Delta, Z, "surv")])
    
    # [T_{(i+1)} - T_{(i)}]
    t_diff <- data_dist[- 1, W] - data_dist[- nrow(data_dist), W]
    
    # Censored subject values (to impute)
    t_cens <- data[data[, Delta] == 0, W]

    # Follow formula assuming Cox model with additional covariates Z
    for (x in which(!uncens)) {
      Zj <- data[x, Z]
      lp <- as.numeric(data.matrix(Zj) %*% matrix(data = fit$coefficients, ncol = 1))
      Cj <- data[x, W]
      Sj <- data_dist[-1, "surv"] ^ (exp(lp)) + data_dist[- nrow(data_dist), "surv"] ^ (exp(lp))
      num <- sum((data_dist[-nrow(data_dist), W] >= Cj) * Sj * t_diff)
      denom <- data[x, "surv"] ^ (exp(lp))
      data$imp[x] <- (1 / 2) * (num / denom) + Cj
    }
  } else {
    if (extrapolate == "efron") {
      # Extend baseline survival curve using Efron's extrapolation 
      extrap_surv0 <- function(t) {
        sapply(X = t, FUN = function(x) {
          if (x < min(surv_df[, W])) {
            1
          } else if (x > max(surv_df[, W])) {
            0
          } else {
            t_before <- which(surv_df[, W] <= x)
            surv_df[max(t_before), "surv0"]
          }
        })
      }
    } else if (extrapolate == "exponential") {
      # Extend baseline survival curve using Brown, Hollander and Kowar's exponential extrapolation 
      extrap_surv0 <- function(t) {
        sapply(X = t, FUN = function(x) {
          if (x < min(surv_df[, W])) {
            1
          } else if (x > max(surv_df[, W])) {
            exp(x * log(surv_df[nrow(surv_df), "surv0"]) / surv_df[nrow(surv_df), W])
          } else {
            t_before <- which(surv_df[, W] <= x)
            surv_df[max(t_before), "surv0"]
          }
        })
      }
    } else if (extrapolate == "weibull") {
      # Estimate Weibull parameters using constrained MLE
      SURVmax <- data[max(which(uncens)), "surv0"]
      Xmax <- data[max(which(uncens)), W] 
      weibull_params <- constr_weibull_mle(t = W, I_event = D, Xtilde = Xmax, rho = SURVmax, alpha0 = 0.1)
      
      # Check that the parameters converged 
      if (!any(is.na(weibull_params))) {
        alpha_hat <- weibull_params[1]
        lambda_hat <- weibull_params[2]
        # Extend baseline survival curve using Moeschberger & Klein's Weibull extrapolation 
        extrap_surv0 <- function(t) {
          sapply(X = t, FUN = function(x) {
            if (x < min(surv_df[, W])) {
              1
            } else if (x > max(surv_df[, W])) {
              exp(- lambda_hat * x ^ alpha_hat)
            } else {
              if (interpolate_between == "carry-forward") {
                # Indices of event times before x
                before <- which(surv_df[, W] < x)
                ## corresponding survival estimate
                surv_df[max(before), "surv0"]
              } else if (interpolate_between == "linear") {
                # Indices of event times before x
                before <- which(surv_df[, W] < x)
                ## Greatest event time before x
                t_before <- surv_df[max(before), W]
                ## corresponding survival estimate
                surv_before <- surv_df[max(before), "surv0"]
                
                # Indices of event times after x
                after <- which(surv_df[, W] > x)
                # Smallest event time after x
                t_after <- surv_df[min(after), W]
                ## corresponding survival estimate
                surv_after <- surv_df[min(after), "surv0"]

                # Linear interpolation of survival estimates before and after
                surv_before + (surv_after - surv_before) / (t_after - t_before) * (x - t_before)
              } else if (interpolate_between == "mean") {
                # Indices of event times before x
                before <- which(surv_df[, W] < x)
                ## corresponding survival estimate
                surv_before <- surv_df[max(before), "surv0"]
                
                # Indices of event times after x
                after <- which(surv_df[, W] > x)
                ## corresponding survival estimate
                surv_after <- surv_df[min(after), "surv0"]
                
                # Mean of survival estimates before and after
                (surv_after + surv_before) / 2
              }
            }
          })
        }
      } else {
        return(list(imputed_data = data, code = FALSE))   
      }
    }  
    
    # Builds on the extrap_surv0 function above by raising S_0(t) ^ HR = S(t|Z)
    extrap_surv <- function(t, hr) {
      sapply(X = t, FUN = extrap_surv0) ^ as.numeric(hr)
    }
    
    ## Use integrate() to approximate integral from W to \infty of S(t|Z)
    int_surv <- sapply(
      X = which(!uncens), 
      FUN = function(i) { 
        tryCatch(expr = integrate(f = extrap_surv, lower = data[i, W], upper = Inf, subdivisions = 2000, hr = data[i, "HR"])$value,
                 error = function(e) return(NA))
      }
    )
    
    # Extrapolate baseline survival S_0(W) for censored W
    ## For \min(X) < W < \max(X), S_0(W) is carried forward from last event 
    ## For W > \max(X), it's extrapolated based on the parameter inputs
    data[which(!uncens), "surv0"] <- sapply(
      X = which(!uncens), 
      FUN = function(i) {
        extrap_surv0(t = data[i, W])
      }
    )
    
    ## And calculate S(t|Z) for censored W 
    data[which(!uncens), "surv"] <- data[which(!uncens), "surv0"] ^ data[which(!uncens), "HR"]
    
    ## Calculate E(X|X>W,Z) = W + int_surv / surv(W|Z)
    data$imp[which(!uncens)] <- data[which(!uncens), W] + int_surv / data[which(!uncens), "surv"]
  }
  
  # Return input dataset with appended column imp containing imputed values 
  if (any(is.na(data$imp))) {
    return(list(imputed_data = data, code = FALSE))  
  } else {
    return(list(imputed_data = data, code = TRUE))
  }
}