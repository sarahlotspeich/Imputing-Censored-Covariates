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

constr_weibull_loglik <- function(alpha, t, I_event, Xtilde, rho) {
  lambda <- - log(rho) / (Xtilde ^ exp(alpha))
  n1 <- sum(I_event)
  weibull_loglik(alpha = alpha, lambda = lambda, t = t, I_event = I_event)
}

constr_weibull_gradient <- function(t, I_event, Xtilde, rho, alpha) {
  n1 <- sum(I_event)
  g <- log(rho) * (Xtilde ^ - alpha) * log(Xtilde) * sum(t ^ alpha)
  g <- g + log(rho) * (Xtilde ^ - alpha) * sum(t ^ alpha * log(t))
  g <- g + sum(I_event * log(t))
  g <- g - n1 * log(Xtilde) + (n1 / alpha)
  return(g)
}

constr_weibull_hessian <- function(t, I_event, Xtilde, rho, alpha) {
  n1 <- sum(I_event)
  H <- log(rho) * (Xtilde ^ (- alpha)) * (log(Xtilde) ^ 2) * sum(t ^ alpha)
  H <- H + log(rho) * (Xtilde ^ (- alpha)) * log(Xtilde) * sum(t ^ alpha * log(t))
  H <- H + log(rho) * (Xtilde ^ (- alpha)) * log(Xtilde) * sum(t ^ alpha * log(t))
  H <- H + log(rho) * (Xtilde ^ (- alpha)) * sum(t ^ alpha * log(t) ^ 2)
  H <- H - n1 / (alpha ^ 2)
  return(H)
}

constr_weibull_mle <- function(t, I_event, Xtilde, rho, alpha0, tol = 1E-4, max_iter = 1E3) {
  it <- 1
  conv <- FALSE
  while (it < max_iter & !conv) {
    g <- constr_weibull_gradient(t = t, I_event = I_event, Xtilde = Xtilde, rho = rho, alpha = alpha0)
    H <- constr_weibull_hessian(t = t, I_event = I_event, Xtilde = Xtilde, rho = rho, alpha = alpha0)
    alpha1 <- alpha0 - g / H
    if (!any(abs(alpha0 - alpha1) > tol)) {
      conv <- TRUE
    }
    it <- it + 1
    alpha0 <- alpha1
  }
  lambda1 <- - log(rho) / (Xtilde ^ alpha1)
  if (conv) {
    return(c(alpha1, lambda1))  
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
    extrap_surv_beyond(x = x, t = t, surv = surv, surv_beyond = beyond, weibull_params = weibull_params)
  }
}