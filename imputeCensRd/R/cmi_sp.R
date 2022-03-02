#' Semiparametric conditional mean imputation (CMI) for a censored predictor
#'
#' Semiparametric conditional mean imputation (CMI) for a censored predictor using a Cox proportional hazards model and the Breslow estimator to estimate conditional survival.
#'
#' @param W Column name of observed predictor values (including censored opens). 
#' @param Delta Column name of censoring indicators. Note that \code{data[, Delta] = 0} is interpreted as a censored observation. 
#' @param Z Column name of additional fully observed covariates.
#' @param data Dataframe or named matrix containing columns \code{W}, \code{Delta}, and \code{Z}.
#' @param fit A \code{coxph} imputation model object modeling \code{W} on \code{Z}. If \code{fit = NULL} (default), the Cox model with only main effects for \code{Z} is fit internally and used.
#' @param extrapolate A string for the method to be used to extrpolate the survival curve beyond the last observed event. Options include \code{"none"} (default), \code{"efron"}, or \code{"brown-hollander-kowar"}.
#'
#' @return A copy of \code{data} with added column \code{imp} containing the imputed values.
#'
#' @export
#' @importFrom survival coxph Surv

cmi_sp <- function(W, Delta, Z, data, fit = NULL, extrapolate = "none") {
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
  surv_df <- with(be, data.frame(t = times, surv = basesurv))
  colnames(surv_df)[1] <- W
  
  # Merge survival estimates into data
  data <- merge(x = data, y = surv_df, all.x = TRUE, sort = FALSE)
  
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
  } else if (extrapolate == "efron") {
    # Extend survival curve using Efron's extrapolation 
    extrap_surv <- function(t) {
      sapply(X = t, FUN = function(x) {
        if (x < min(surv_df[, W])) {
          1
        } else if (x > max(surv_df[, W])) {
          0
        } else {
          t_before <- which(surv_df[, W] <= x)
          surv_df[max(t_before), "surv"]
        }
      })
    }
    ## Use integrate() to approximate integral from W to \infty of S(t|Z)
    int_surv <- sapply(
      X = which(!uncens), 
      FUN = function(i) { 
        tryCatch(expr = integrate(f = extrap_surv, lower = data[i, W], upper = Inf, subdivisions = 2000)$value,
                 error = function(e) return(NA))
      }
    )
  } else if (extrapolate == "gill") {
    # Extend survival curve using Gill's extrapolation 
    extrap_surv <- function(t) {
      sapply(X = t, FUN = function(x) {
        if (x < min(surv_df[, W])) {
          1
        } else if (x > max(surv_df[, W])) {
          surv_df[nrow(surv_df), "surv"]
        } else {
          t_before <- which(surv_df[, W] <= x)
          surv_df[max(t_before), "surv"]
        }
      })
    }
  } else if (extrapolate == "brown-hollander-kowar") {
    # Extend survival curve using Gill's extrapolation 
    extrap_surv <- function(t, hr) {
      sapply(X = t, FUN = function(x, hr) {
        if (x < min(surv_df[, W])) {
          1
        } else if (x > max(surv_df[, W])) {
          exp(x * log(surv_df[nrow(surv_df), "surv"] ^ hr) / surv_df[nrow(surv_df), W])
        } else {
          t_before <- which(surv_df[, W] <= x)
          surv_df[max(t_before), "surv"] ^ hr
        }
      }, hr = hr)
    }
    ## Use integrate() to approximate integral from W to \infty of S(t|Z)
    int_surv <- sapply(
      X = which(!uncens), 
      FUN = function(i) { 
        tryCatch(expr = integrate(f = extrap_surv, lower = data[i, W], upper = Inf, subdivisions = 2000, hr = data[i, "HR"])$value,
                 error = function(e) return(NA))
      }
    )
  }
  
  if (extrapolate %in% c("efron", "gill", "brown-hollander-kowar")) {
    ## Calculate E(X|X>W,Z) = int_surv / surv(W|Z) + W
    data$imp[which(!uncens)] <- data[which(!uncens), W] + int_surv / data[which(!uncens), "surv"]
  }
  
  # Return input dataset with appended column imp containing imputed values 
  return(data)
}