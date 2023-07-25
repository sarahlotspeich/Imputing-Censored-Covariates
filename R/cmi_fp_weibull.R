#' @export
cmi_fp_weibull = function(W, Delta, Z, data, fit, infinite_integral = TRUE, maxiter = 100) {
  # Initialize imputed values 
  data$imp = data[, W] ## start with imp = W
  
  # Calculate linear predictor \lambda_0 + \blambda_1 %*% Z for AFT model
  lp = fit$coefficients[1] +
    data.matrix(data[, Z]) %*% matrix(data = fit$coefficients[- 1], ncol = 1) ## linear predictors
  
  # Transform parameters to agree with R's weibull parameterization 
  weib_shape = 1 / fit$scale
  weib_scale = exp(lp)
  
  # Create an indicator variable for being uncensored
  uncens = data[, Delta] == 1
  
  if (infinite_integral) {
    # Use adaptive quadrature to estimate
    ## integral from X = 0 to X = Wi
    int_surv_W_to_Inf = sapply(
      X = which(!uncens), 
      FUN = function(i) { 
        integrate(f = pweibull, 
                  lower = data[i, W], 
                  upper = Inf, 
                  shape = weib_shape, 
                  scale = weib_scale[i],
                  lower.tail = FALSE)$value
      }
    )

    
  } else {
    # Calculate mean life = integral from 0 to \infty of S(t|Z)
    est_ml = weib_scale[which(!uncens)] * gamma(1 + fit$scale) ## using formula for Weibull distribution
    
    # Use adaptive quadrature to estimate
    ## integral from X = 0 to X = Wi
    int_surv_0_to_W = sapply(
      X = which(!uncens), 
      FUN = function(i) { 
        integrate(f = pweibull, 
                  lower = 0, 
                  upper = data[i, W], 
                  shape = weib_shape, 
                  scale = weib_scale[i],
                  lower.tail = FALSE)$value
      }
    )
    
    # Subtract integral from mean life to get 
    ## integral from X = Wi to X = Infinity
    int_surv_W_to_Inf = est_ml - int_surv_0_to_W
  }
  
  ## Calculate S(W|Z)
  surv = pweibull(q = data[which(!uncens), W], 
                  shape = weib_shape, 
                  scale = weib_scale[which(!uncens)], 
                  lower.tail = FALSE)
  
  ## Calculate MRL(W) = int_surv / S(W|Z)
  est_mrl = int_surv_W_to_Inf / surv
  
  ## Calculate imputed value E(X|X>W,Z) = W + MRL(W)
  data[which(!uncens), "imp"] = data[which(!uncens), W] + est_mrl
  
  ## Check for infinite/NA imputed values 
  # if (any(is.na(data$imp))) {
  #   data$imp[which(is.na(data$imp))] = data[which(is.na(data$imp)), W]
  # }
  # if (any(data$imp == Inf)) {
  #   data$imp[which(data$imp == Inf)] = data[which(data$imp == Inf), W]
  # }
  
  # Return input dataset with appended column imp containing imputed values 
  return(list(imputed_data = data, code = !any(is.na(data$imp))))  
}