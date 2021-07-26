# n: sample size for single simulation
# n.sims: number of simulations
# beta0: true intercept
# beta1: true slope on x
# sigma: true sd for epsilon
# xshape: shape parameter for sim. of x
# xscale: scale parameter for sim. of x
# cshape: shape parameter for sim. of c
# cscale: scale parameter for sim. of c
generate.data = function(n, n.sims, 
                         beta0, betaX, betaZ = NULL, sigma = 1,
                         x = NULL, z = NULL,
                         xshape = 0.75, xscale = 0.25,
                         c.lower = 0, c.upper = 1) {
  # total number of rows
  N = n*n.sims
  
  # subject id and simulation id
  subj.id <- rep(1:n, n.sims)
  sim.id <- rep(1:n.sims, each = n)

  # Data generation
  # Generate x if not provided
  if (is.null(x)) { x <- rweibull(n = N, shape = xshape, scale = xscale) }
  
  # Generate random error and outcome
  epsilon <- rnorm(n = N, sd = sigma)
  y <- beta0 + betaX * x  + epsilon
  
  # Censoring mechanism
  cens <- runif(n = N, min = c.lower, max = c.upper)
  delta <- as.numeric(x <= cens)
  t <- ifelse(delta, x, cens)
  
  # If no Z variables, return id's, response, t, and delta
  if (is.null(betaZ)) {
    data.frame(subj.id, sim.id, y, t, delta) %>% 
      return()
  }
  else {
    pZ <- length(betaZ)
    # Generate z if not provided
    if (is.null(z)) { z <- rnorm(n = N*pZ, mean = 0, sd = 1) %>% matrix(nrow = N, ncol = pZ) }
    colnames(z) = paste0("Z", 1:pZ)
    y <- y + z %*% betaZ
    data.frame(subj.id, sim.id, y, t, delta) %>% 
      cbind(z) %>%
      return()
  }
}

# Simulate survival data based on Bender (2005)
cox.simulation = function(n, param, covariate, dist = "Exponential", 
                          lambda = 1, nu = NULL, alpha = NULL) {
  U = runif(n = n, min = 0, max = 1)
  lin.pred = covariate %*% param
  H_0.T = -log(U)*exp(-lin.pred)
  
  if (dist == "Exponential") { return(H_0.T/lambda) }
  else if (dist == "Weibull") { return((H_0.T/lambda)^{1/nu}) }
  else if (dist == "Gompertz") { return(log(1 + alpha*t/lambda)/alpha) }
  else { print("Invalid distribution") }
}

# Generates accelerated failure time (AFT) response
# n: desired simulation sample size
# gamma: parameter(s) for aft model
# a: matrix of covariates for aft model, genertaed from N(0, 1) if not provided
# sigma: standard deviation of random error
aft.simulation = function(n, gamma, a = NULL, sigma = 0.1) {
  # number of parameters in AFT model
  p <- length(gamma)
  
  # generate auxiliary variables if not provided
  if (is.null(a)) { 
    a <- rnorm(n = n*(p - 1), sd = 1) %>%
      matrix(nrow = n, ncol = p - 1)
  }
  
  # generate outcome
  x <- exp(gamma[1] + a %*% gamma[-1] + rnorm(n = n, sd = sigma))
  
  # data table to be returned as result
  result = data.table(x = x)
  result = result %>%
    cbind(as.data.table(a))
  
  # name columns of data table and return
  colnames(result) = c("x", paste0("A",1:(p - 1)))
  return(result)
}

# Returns vector of preceding events (x[i] for which delta[i] == 1)
# for each element of x
# x: vector of event times
# delta: vector of event indicators
prev.event = function(x, delta) {
  result = x
  events = x[delta == 1]
  for (i in which(delta == 0)) {
    result[i] = max(events[events < result[i]])
  }
  
  return(result)
}