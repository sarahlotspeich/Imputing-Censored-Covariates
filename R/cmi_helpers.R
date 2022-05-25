surv_below_minX <- function(data, W, Delta, S) {
  # Create an indicator variable for being uncensored
  uncens <- data[, Delta] == 1
  
  # Assume survival at censored W < min(X) = 1
  minX <- as.numeric(data[min(which(uncens)), W]) 
  data[which(data[, W] < minX), S] <- 1
  # Return data with filled in survival
  return(data)
}

surv_between_X <- function(data, W, Delta, S, surv_between) {
  # Create an indicator variable for being uncensored
  uncens <- data[, Delta] == 1
  
  # Save data frame with survival estimates (ordered)
  surv_df <- data[data[, Delta] == 1, ]
  
  # Interpolate survival at censored W < \widetilde{X}
  Xtilde <- as.numeric(data[max(which(uncens)), W])
  needs_interp <- which(is.na(data[, S]) & data[, W] < Xtilde)
  data[needs_interp, S] <- sapply(X = data[needs_interp, W], 
                                  FUN = interp_surv_between, 
                                  t = surv_df[, W], 
                                  surv = surv_df[, S], 
                                  surv_between = surv_between)
  
  # Return data with interpolated survival
  return(data)
}

surv_beyond_maxX <- function(data, W, Delta, S, Xmax, surv_beyond) {
  # Create an indicator variable for being uncensored
  uncens <- data[, Delta] == 1
  
  # Save data frame with survival estimates (ordered)
  surv_df <- data[data[, Delta] == 1, ]
  
  # Where to tie into the survival estimator
  Xtilde <- max(data[uncens, W]) 
  rho <- min(data[uncens, S], na.rm = TRUE) 
  
  # Weibull parameters 
  if (surv_beyond == "weibull") {
    if (Xmax < Inf) {
      weibull_params <- dbl_constr_weibull(Xtilde = Xtilde, 
                                           rho = rho, 
                                           Xmax = Xmax)
    } else {
      weibull_params <- constr_weibull_mle(t = data[, W], 
                                           I_event = data[, Delta], 
                                           Xtilde = Xtilde, 
                                           rho = rho, 
                                           alpha0 = 1E-4)
    }
    
    # If weibull params don't converge, quit 
    if (any(is.na(weibull_params))) {
      return(list(data = data, code = FALSE, weibull_params = c(NA, NA)))   
    }
  } else {
    weibull_params <- NULL
  }
  
  # Extrapolate survival at censored W > \widetilde{X}
  needs_extrap <- which(is.na(data[, S]) & data[, W] > Xtilde)
  data[needs_extrap, S] <- sapply(X = data[needs_extrap, W], 
                                  FUN = extrap_surv_beyond, 
                                  t = surv_df[, W], 
                                  surv = surv_df[, S], 
                                  surv_beyond = surv_beyond, 
                                  weibull_params = weibull_params)
  
  # Return data with extrapolated survival (+ convergence code for Weibull)
  return(list(data = data, code = TRUE, weibull_params = weibull_params))   
}

impute_cens_W <- function(data, W, Delta, Z, S, S0, Xmax, trapezoidal_rule, surv_between, surv_beyond, weibull_params) {
  # Save data frame with survival estimates (ordered)
  surv_df <- data[data[, Delta] == 1, ]
  
  # For people with events, imp = X
  data$imp <- data[, W]
  
  # Row IDs in need of imputation
  needs_impute <- which(data[, Delta] == 0) 
  
  # Calculate imputed values 
  if (trapezoidal_rule) {
    # Distinct rows (in case of non-unique obs values)
    data_dist <- unique(data[, c(W, Delta, Z, S)])
    
    # [T_{(i+1)} - T_{(i)}]
    t_diff <- data_dist[- 1, W] - data_dist[- nrow(data_dist), W]
    
    # S(t+1) + S(t)
    surv_sum <- data_dist[-1, S] + data_dist[- nrow(data_dist), S]
    
    # Use trapezoidal approximation for integral
    for (i in needs_impute) {
      sum_surv_i <- sum((data_dist[- nrow(data_dist), W] >= as.numeric(data[i, W])) * surv_sum * t_diff)
      data$imp[i] <- data$imp[i] + (1 / 2) * (sum_surv_i /  data[i, S])
    }
  } else {
    
    if (!is.null(S0)) {
      # Builds on the extend_surv function by raising S_0(t) ^ HR = S(t|Z)
      to_integrate <- function(t, hr) {
        basesurv <- sapply(X = t, 
                           FUN = extend_surv, 
                           t = surv_df[, W], 
                           surv = surv_df[, S0], 
                           surv_between = surv_between, 
                           surv_beyond = surv_beyond, 
                           weibull_params = weibull_params)  
        basesurv ^ as.numeric(hr)
      }
      
      ## Use integrate() to approximate integral from W to Xmax of S(t|Z)
      int_surv <- sapply(
        X = needs_impute, 
        FUN = function(i) { 
          tryCatch(expr = integrate(f = to_integrate, 
                                    lower = data[i, W], 
                                    upper = Xmax, 
                                    subdivisions = 2000, 
                                    hr = data[i, "HR"])$value,
                   error = function(e) return(NA))
        }
      )
    } else {
      to_integrate <- function(t) {
        surv <- sapply(X = t, 
                       FUN = extend_surv, 
                       t = surv_df[, W], 
                       surv = surv_df[, S],
                       surv_between = surv_between, 
                       surv_beyond = surv_beyond, 
                       weibull_params = weibull_params)  
        surv
      } 
      
      ## Use integrate() to approximate integral from W to Xmax of S(t)
      int_surv <- sapply(
        X = needs_impute, 
        FUN = function(i) { 
          tryCatch(expr = integrate(f = to_integrate, 
                                    lower = data[i, W], 
                                    upper = Xmax, 
                                    subdivisions = 2000)$value,
                   error = function(e) return(NA))
        }
      )
    }
    
    ## Calculate E(X|X>W) = W + int_surv / surv(W)
    data$imp[needs_impute] <- data[needs_impute, W] + int_surv / data[needs_impute, S]
  }
  
  # Return data with imputed predictors (imp)
  return(data)
}

extend_and_impute <- function(data, W, Delta, Z, S, S0 = NULL, Xmax, trapezoidal_rule, surv_between, surv_beyond) {
  # Create an indicator variable for being uncensored
  uncens <- data[, Delta] == 1
  
  # Assume survival at censored W < min(X) = 1
  # minX <- data[min(which(uncens & has_z)), W] 
  # data[which(data[, W] < minX & has_z), "surv"] <- 1
  data <- surv_below_minX(data = data, 
                          W = W, 
                          Delta = Delta, 
                          S = S)
  
  # Interpolate survival at censored W < \widetilde{X}
  # Xtilde <- data[max(which(uncens & has_z)), W] 
  # needs_interp <- which(is.na(data[, "surv"]) & has_z & data[, W] < Xtilde)
  # data[needs_interp, "surv"] <- sapply(X = data[needs_interp, W], 
  #                                      FUN = interp_surv_between, 
  #                                      t = km_surv[which(km_surv[, Z] == z), W], 
  #                                      surv = km_surv[which(km_surv[, Z] == z), "surv"], 
  #                                      surv_between = surv_between)
  data <- surv_between_X(data = data, 
                         W = W, 
                         Delta = Delta, 
                         S = S, 
                         surv_between = surv_between)
  
  # Save parameters to tie Weibull curve into KM estimator
  Xtilde <- max(data[uncens, W]) 
  rho <- min(data[uncens, S], na.rm = TRUE) 
  
  # Extrapolate survival at censored W > \widetilde{X}
  if (rho > 0) {
    extrap <- surv_beyond_maxX(data = data, 
                               W = W, 
                               Delta = Delta, 
                               S = S, 
                               Xmax = Xmax, 
                               surv_beyond = surv_beyond)

    # (If surv_beyond = "Weibull")
    ## Return without imputing if parameters did not converge
    if (!extrap$code) {
      data$imp <- NA
      return(list(data = data, code = FALSE))   
    }
    
    # Calculate imputed values E(X|X>W,Z)
    data <- impute_cens_W(data = extrap$data, 
                          W = W, 
                          Delta = Delta, 
                          Z = Z, 
                          S = S, 
                          S0 = S0, 
                          Xmax = Xmax, 
                          trapezoidal_rule = trapezoidal_rule, 
                          surv_between = surv_between, 
                          surv_beyond = surv_beyond, 
                          weibull_params = extrap$weibull_params)
  } else {
    # If survival estimator ends at 0, no need to extrapolate
    ## For computational ease, treat this as "immediate drop off"
    data <- impute_cens_W(data = data, 
                          W = W, 
                          Delta = Delta, 
                          Z = Z, 
                          S = S, 
                          S0 = S0, 
                          Xmax = Xmax, 
                          trapezoidal_rule = trapezoidal_rule, 
                          surv_between = surv_between, 
                          surv_beyond = "drop", 
                          weibull_params = NA)
  }
  
  return(list(data = data, code = TRUE))
}