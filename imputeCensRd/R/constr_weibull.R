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
