solve_for_shape <- function(scale, SURVmax, Xmax) {
  find_root <- uniroot(f = function(shape) calc_weibull_surv(t = Xmax, scale = scale, shape = shape) - SURVmax, 
                       interval = c(0.1, 2), extendInt = "yes")
  return(find_root$root)
}

calc_weibull_surv <- function(t, scale, shape) {
  exp(- (t / scale) ^ shape)
}

constr_weib_loglik <- function(scale, W, Delta, data) {
  ####################################################
  # Pre-processing ###################################
  ####################################################
  # < number of uncensored subjects > ----------------
  n1 <- sum(data[, Delta]) # -------------------------
  # ---------------- < number of uncensored subjects >
  # Reordered data by W ------------------------------
  data <- data[order(data[, W], decreasing = FALSE), ]
  # ------------------------------ Reordered data by W
  # Reordered data to be uncensored first ------------
  data <- data[order(data[, Delta], decreasing = TRUE), ]
  # ------------ Reordered data to be uncensored first
  # Create subset of uncensored subjects' data -------
  uncens_data <- data[1:n1, ]
  # ------- Create subset of uncensored subjects' data
  # Create subset of censored subjects' data ---------
  cens_data <- data[-c(1:n1), ]
  # --------- Create subset of censored subjects' data
  
  # Use constraint to find other parameter -----------
  shape <- solve_for_shape(scale = scale, 
                           Xmax = uncens_data[nrow(uncens_data), W], 
                           SURVmax = uncens_data[nrow(uncens_data), "surv0"])
  # ----------- Use constraint to find other parameter
  
  # Calculate survival at uncensored -----------------
  sXgivZ <- calc_weibull_surv(scale = scale, shape = shape, t = uncens_data[, W])
  # ----------------- Calculate survival at uncensored
  
  ####################################################
  # Calculate the log-likelihood #####################
  ####################################################
  # Log-likelihood contribution of uncensored X ------
  ll <- sum(log(sXgivZ))
  # ------ Log-likelihood contribution of uncensored X
  # Log-likelihood contribution of censored X --------
  integrate_surv <- function(data_row) {
    data_row <- data.frame(t(data_row))
    return(
      tryCatch(expr = integrate(f = calc_weibull_surv, lower = data_row[, W], upper = Inf, 
                                scale = scale, shape = shape)$value,
               error = function(err) {0})
    )
  }
  integral <- apply(X = cens_data, MARGIN = 1, FUN = integrate_surv)
  log_integral <- log(integral)
  log_integral[log_integral == -Inf] <- 0
  ll <- ll + sum(log_integral)
  # -------- Log-likelihood contribution of censored X
  # Return (-1) x log-likelihood for use with nlm() --
  return(- ll)
  # -- Return (-1) x log-likelihood for use with nlm()
}

get_weib_params <- function(W, Delta, data) {
  ####################################################
  # Pre-processing ###################################
  ####################################################
  # < number of uncensored subjects > ----------------
  n1 <- sum(data[, Delta]) # -------------------------
  # ---------------- < number of uncensored subjects >
  # Reordered data by W ------------------------------
  data <- data[order(data[, W], decreasing = FALSE), ]
  # ------------------------------ Reordered data by W
  # Reordered data to be uncensored first ------------
  data <- data[order(data[, Delta], decreasing = TRUE), ]
  # ------------ Reordered data to be uncensored first
  # Create subset of uncensored subjects' data -------
  uncens_data <- data[1:n1, ]
  # ------- Create subset of uncensored subjects' data
  # Create subset of censored subjects' data ---------
  cens_data <- data[-c(1:n1), ]
  # --------- Create subset of censored subjects' data
  
  # try_scale <- seq(0.1, 10, by = 0.1)
  # Use constraint to find other parameter -----------
  # try_shape <- sapply(X = try_scale, FUN = solve_for_shape, 
  #                    Xmax = uncens_data[nrow(uncens_data), W], 
  #                    SURVmax = uncens_data[nrow(uncens_data), "surv0"])
  # ----------- Use constraint to find other parameter
  
  suppressWarnings(
    fit <- nlm(f = constr_weib_loglik, p = 0.1, W = W, Delta = Delta, data = data)
  )
  
  if (fit$code <= 2) {
    fit_scale <- fit$estimate
    n1 <- sum(data[, Delta])
    fit_shape <- solve_for_shape(scale = fit_scale, 
                                 Xmax = uncens_data[nrow(uncens_data), W], 
                                 SURVmax = uncens_data[nrow(uncens_data), "surv0"])
    return(c(fit_scale, fit_shape))  
  }
}
