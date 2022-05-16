interp_surv_between <- function(x, t, surv, surv_between) {
  if (surv_between == "carry-forward") {
    # Indices of event times before x
    before <- which(t <= x)
    ## corresponding survival estimate
    surv[max(before)]
  } else if (surv_between == "linear") {
    # Indices of event times before x
    before <- which(t <= x)
    ## Greatest event time before x
    t_before <- t[max(before)]
    ## corresponding survival estimate
    surv_before <- surv[max(before)]
    # Indices of event times after x
    after <- which(t >= x)
    # Smallest event time after x
    t_after <- t[min(after)]
    ## corresponding survival estimate
    surv_after <- surv[min(after)]
    # Linear interpolation of survival estimates before and after
    surv_before + (surv_after - surv_before) / (t_after - t_before) * (x - t_before)
  } else if (surv_between == "mean") {
    # Indices of event times before x
    before <- which(t <= x)
    ## corresponding survival estimate
    surv_before <- surv[max(before)]
    # Indices of event times after x
    after <- which(t >= x)
    ## corresponding survival estimate
    surv_after <- surv[min(after)]
    # Mean of survival estimates before and after
    (surv_after + surv_before) / 2
  }
}

weibull_loglik <- function(alpha, lambda, t, I_event) {
  n1 <- sum(I_event)
  ll <- - lambda * sum(t ^ alpha)
  ll <- ll + (alpha - 1) * sum(I_event * log(t)) 
  ll <- ll + n1 * log(lambda)
  ll <- ll + n1 * log(alpha)
  return(- ll)
}

constr_weibull_loglik <- function(alpha, t, I_event, Xtilde, rho) {
  lambda <- - log(rho) / (Xtilde ^ exp(alpha))
  n1 <- sum(I_event)
  weibull_loglik(alpha = alpha, lambda = lambda, t = t, I_event = I_event)
}

constr_weibull_mle <- function(t, I_event, Xtilde, rho, alpha0, tol = 1E-4, max_iter = 1E3) {
  t[t == 0] <- 1E-4
  suppressWarnings(
    nlm_res <- nlm(f = constr_weibull_loglik, 
                   p = alpha0, 
                   t = t, 
                   I_event = I_event, 
                   Xtilde = Xtilde, 
                   rho = rho)
  )
  conv <- nlm_res$code <= 2 & nlm_res$iterations > 0 
  if (conv) {
    alpha1 <- nlm_res$estimate
    lambda1 <- - log(rho) / (Xtilde ^ alpha1)
    return(c(alpha1, lambda1))  
  } else {
    return(c(NA, NA))
  }
}

surv_omega <- function(alpha, Xtilde, rho, Xmax, approx0) {
  lambda <- - log(rho) / Xtilde ^ alpha
  exp(- lambda * Xmax ^ alpha) - approx0
}

dbl_constr_weibull <- function(Xtilde, rho, Xmax, alpha_l = 1E-4, alpha_u = 10, tol = 1E-2, approx0 = 1E-4) {
  find_root <- tryCatch(expr = uniroot(f = surv_omega, 
                                       lower = alpha_l, 
                                       upper = alpha_u, 
                                       extendInt = "yes",
                                       rho = rho,
                                       Xtilde = Xtilde, 
                                       Xmax = Xmax,
                                       approx0 = approx0),
                        error = function(e) return(list(estim.prec = 99999)))
  if (find_root$estim.prec < tol) {
    alpha_c <- find_root$root
    lambda_c <- - log(rho) / Xtilde ^ alpha_c
    return(c(alpha_c, lambda_c))
  } else {
    return(c(NA, NA))
  }
}

extrap_surv_beyond <- function(x, t, surv, surv_beyond, weibull_params = NULL) {
  if (surv_beyond == "carry-forward") {
    before <- which(t <= x)
    surv[max(before)]
  } else if (surv_beyond == "drop-off") {
    0 
  } else if (surv_beyond == "exponential") {
    exp(x * log(surv[length(surv)]) / t[length(t)])
  } else if (surv_beyond == "weibull") {
    alpha_hat <- weibull_params[1]
    lambda_hat <- weibull_params[2]
    exp(- lambda_hat * x ^ alpha_hat)
  }
}

extend_surv <- function(x, t, surv, surv_between, surv_beyond, weibull_params = NULL) {
  Xtilde <- max(t)
  if (x < min(t)) {
    ## if x <= min(t), 100% survival
    1
  } else if (x <= Xtilde) {
    ## if x <= Xtilde, use interpolate functions
    interp_surv_between(x = x, t = t, surv = surv, surv_between = surv_between)
  } else {
    ## if x > Xtilde, use extrapolate functions
    extrap_surv_beyond(x = x, t = t, surv = surv, surv_beyond = surv_beyond, weibull_params = weibull_params)
  }
}