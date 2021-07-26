# n: sample size for single simulation
# n.sims: number of simulations
# beta0: intercept
# betaX: coefficient for covariate x
# betaZ: coefficient for covariate(s) z
# sigma: sd for random error in outcome model
# xshape: shape parameter for sim. of x
# xscale: scale parameter for sim. of x
# c.lower: minimum for sim. of c
# c.upper: maximum for sim. of c
generate.data = function(n, n.sims, 
                         beta0, betaX, betaZ = NULL, sigma = 1,
                         x = NULL, z = NULL,
                         xshape = 0.75, xscale = 0.25,
                         c.lower = 0, c.upper = 1) {
  # total number of rows
  N <- n*n.sims
  
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
  
  # If no z coefficient is provided, return id's, response, t, and delta
  if (is.null(betaZ)) {
    data.frame(subj.id, sim.id, y, t, delta) %>% 
      return()
  }
  else {
    pZ <- length(betaZ)
    # Generate design matrix fof covariates z if not provided
    if (is.null(z)) { z <- rnorm(n = N*pZ, mean = 0, sd = 1) %>% matrix(nrow = N, ncol = pZ) }
    colnames(z) = paste0("Z", 1:pZ)
    # Add the effect of z to the outcome
    y <- y + z %*% betaZ
    # return id's, response, t, delta, and z
    data.frame(subj.id, sim.id, y, t, delta) %>% 
      cbind(z) %>%
      return()
  }
}

# Simulate survival data based on Bender (2005)
# n: sample size
# param: coefficients for simulation
# covariate: covaraites for simulation
# dist: desired outcome distribution
# lambda: simulation parameter
# nu: simulation parameter
# alpha: simulation parameter
cox.simulation = function(n, param, covariate, dist = "Exponential", 
                          lambda = 1, nu = NULL, alpha = NULL) {
  # Generate n iid unif(0, 1) random variable
  U = runif(n = n, min = 0, max = 1)
  # Calculate H_0(T) according to Bender (2005)
  H_0.T = -log(U)*exp(-covariate %*% param)
  
  # Return simulated response based on distribution character provided
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
  
  # generate covariates if not provided
  if (is.null(a)) { 
    a <- rnorm(n = n*(p - 1), sd = 1) %>%
      matrix(nrow = n, ncol = p - 1)
  }
  
  # generate outcome
  x <- exp(gamma[1] + a %*% gamma[-1] + rnorm(n = n, sd = sigma))
  
  # data frame to be returned as result
  result = data.frame(x = x)
  result = result %>%
    cbind(as.data.frame(a))
  
  # name columns of data and return
  colnames(result) = c("x", paste0("A",1:(p - 1)))
  return(result)
}