#' @export
cmi_fp_weibull = function(W, Delta, Z, data, fit, infinite_integral = TRUE, maxiter = 100) {
  # Initialize imputed values 
  data$imp = data[, W] ## start with imp = W
  
  # Calculate linear predictor \lambda_0 + \blambda_1 %*% Z for AFT model
  lp = fit$coefficients[1] +
    data.matrix(data[, Z]) %*% matrix(data = fit$coefficients[- 1], ncol = 1) ## linear predictors
  
  # Create an indicator variable for being uncensored
  uncens = data[, Delta] == 1
  
  if (infinite_integral) {
    # Use adaptive quadrature to estimate
    ## integral from X = 0 to X = Wi
    int_surv_W_to_Inf = sapply(
      X = which(!uncens), 
      FUN = function(i) { 
        integrate(f = my_psurvreg, 
                  lower = data[i, W], 
                  upper = Inf, 
                  mean = lp[i], 
                  scale = fit$scale,
                  distribution = "weibull",
                  lower.tail = FALSE)$value
      }
    )
    
    ## Calculate S(W|Z)
    surv = my_psurvreg(q = data[which(!uncens), W],
                       mean = lp[which(!uncens)], 
                       scale = fit$scale, 
                       distribution = "weibull", 
                       lower.tail = FALSE)
    
    ## Calculate MRL(W) = int_surv / S(W|Z)
    est_mrl = int_surv_W_to_Inf / surv
  } else {
    # Calculate mean life = integral from 0 to \infty of S(t|Z)
    est_ml = exp(lp[which(!uncens)]) * gamma(1 + fit$scale) ## using formula for Weibull distribution
    
    # Use adaptive quadrature to estimate
    ## integral from X = 0 to X = Wi
    int_surv_0_to_W = sapply(
      X = which(!uncens), 
      FUN = function(i) { 
        integrate(f = my_psurvreg, 
                  lower = 0, 
                  upper = data[i, W], 
                  mean = lp[i], 
                  scale = fit$scale,
                  distribution = "weibull",
                  lower.tail = FALSE)$value
      }
    )
    
    ## Calculate imputed value 
    int_surv_W_to_Inf = est_ml - int_surv_0_to_W
    surv = my_psurvreg(q = data[which(!uncens), W],
                       mean = lp[which(!uncens)], 
                       scale = fit$scale, 
                       distribution = "weibull", 
                       lower.tail = FALSE)
    
    ## Calculate MRL(W) = int_surv / S(W|Z)
    est_mrl = int_surv_W_to_Inf / surv
  }
  
  ## Calculate imputed value E(X|X>W,Z) = W + MRL(W)
  data[which(!uncens), "imp"] = data[which(!uncens), W] + est_mrl
  
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