
impute_cens_W_km <- function(data, W, Delta, Z, S, Xmax, trapezoidal_rule, surv_between, surv_beyond, weibull_params) {
  # Save data frame with survival estimates (ordered)
  surv_df <- data[data[, Delta] == 1, ]
  
  # For people with events, imp = X
  data$imp <- data[, W]
  
  # Row IDs in need of imputation
  needs_impute <- which(data[, Delta] == 0) 
  
  # Check that imputation is needed
  if (length(needs_impute) > 0) {
    # Calculate imputed values 
    if (trapezoidal_rule) {
      # Distinct rows (in case of non-unique obs values)
      data_dist <- unique(data[, c(W, Delta, Z, S)])
      
      # [T_{(i+1)} - T_{(i)}]
      t_diff <- data_dist[- 1, W] - data_dist[- nrow(data_dist), W]
      
      # S(t+1) + S(t)
      surv_sum <- data_dist[- 1, S] + data_dist[- nrow(data_dist), S]
      
      # Use trapezoidal approximation for integral
      for (i in needs_impute) {
        sum_surv_i <- sum((data_dist[- nrow(data_dist), W] >= as.numeric(data[i, W])) * surv_sum * t_diff)
        data$imp[i] <- data$imp[i] + (1 / 2) * (sum_surv_i /  data[i, S])
      }
    } else {
      # Integrand (survival function)
      to_integrate <- function(t) {
        sapply(X = t, 
               FUN = extend_surv, 
               t = surv_df[, W], 
               surv = surv_df[, S],
               surv_between = surv_between, 
               surv_beyond = surv_beyond, 
               weibull_params = weibull_params) 
      } 
      
      # Use integrate() to approximate integral from W to Xmax of S(t)
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

extend_and_impute_km <- function(data, W, Delta, Z, S, Xmax, trapezoidal_rule, surv_between, surv_beyond) {
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
    data <- impute_cens_W_km(data = extrap$data, 
                             W = W, 
                             Delta = Delta, 
                             Z = Z, 
                             S = S, 
                             Xmax = Xmax, 
                             trapezoidal_rule = trapezoidal_rule, 
                             surv_between = surv_between, 
                             surv_beyond = surv_beyond, 
                             weibull_params = extrap$weibull_params)
  } else {
    # If survival estimator ends at 0, no need to extrapolate
    ## For computational ease, treat this as "immediate drop off"
    data <- impute_cens_W_km(data = data, 
                             W = W, 
                             Delta = Delta, 
                             Z = Z, 
                             S = S,
                             Xmax = Xmax, 
                             trapezoidal_rule = trapezoidal_rule, 
                             surv_between = surv_between, 
                             surv_beyond = "drop", 
                             weibull_params = NA)
  }
  
  return(list(data = data, code = TRUE))
}
