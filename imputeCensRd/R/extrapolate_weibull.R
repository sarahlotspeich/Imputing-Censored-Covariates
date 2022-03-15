solve_for_shape <- function(scale, SURVmax, Xmax) {
  find_root <- uniroot(f = function(shape) calc_weibull_surv(t = Xmax, scale = scale, shape = shape) - SURVmax, 
                       interval = c(0.1, 2), extendInt = "yes")
  return(find_root$root)
}

solve_for_scale <- function(shape, SURVmax, Xmax) {
  scale <- - Xmax ^ shape / log(SURVmax)
  return(scale)
}

calc_weibull_surv <- function(t, scale, shape) {
  exp(- (t / scale) ^ shape)
}

calc_weibull_surv_k <- function(t, k, Xmax, SURVmax) {
  exp((t / Xmax) ^ k * log(SURVmax))
}

constr_weib_loglik <- function(k, W, Delta, data) {
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
  
  SURVmax <- uncens_data[nrow(uncens_data), "surv0"]
  Xmax <- uncens_data[nrow(uncens_data), W]
  
  # Calculate survival at uncensored -----------------
  sXgivZ <- calc_weibull_surv_k(t = uncens_data[, W], k = k, Xmax = Xmax, SURVmax = SURVmax)
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
      tryCatch(expr = integrate(f = calc_weibull_surv_k, lower = data_row[, W], upper = Inf, 
                                Xmax = Xmax, SURVmax = SURVmax)$value,
               error = function(err) {0})
    )
  }
  integral <- apply(X = cens_data, MARGIN = 1, FUN = integrate_surv)
  log_integral <- log(integral)
  log_integral[log_integral == -Inf] <- 0
  ll <- ll + sum(log_integral)
  # -------- Log-likelihood contribution of censored X
  # Return log-likelihood for use with nlm() ---------
  return(- ll)
  # ------- Return x log-likelihood for use with nlm()
}

get_weib_params <- function(W, Delta, data) {
  suppressWarnings(
    fit <- nlm(f = constr_weib_loglik, p = 0.1, W = W, Delta = Delta, data = data)
  )
  
  if (fit$code <= 2) {
    return(fit$estimate)  
  }  
}