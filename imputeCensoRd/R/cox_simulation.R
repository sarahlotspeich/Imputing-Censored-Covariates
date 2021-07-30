# Simulate survival data based on Bender (2005)
# n: sample size
# logHR: coefficients for simulation
# covariate: covaraites for simulation
# dist: desired outcome distribution
# lambda: simulation parameter
# nu: simulation parameter
# alpha: simulation parameter
cox.simulation = function(n, logHR, covariate, dist = "Exponential", 
                          lambda = 1, nu = NULL, alpha = NULL) {
  # Generate n iid unif(0, 1) random variable
  U = runif(n = n, min = 0, max = 1)
  # Calculate H_0(T) according to Bender (2005)
  H_0.T = -log(U)*exp(-covariate %*% logHR)
  
  # Return simulated response based on distribution character provided
  if (dist == "Exponential") { return(H_0.T/lambda) }
  else if (dist == "Weibull") { return((H_0.T/lambda)^{1/nu}) }
  else if (dist == "Gompertz") { return(log(1 + alpha*t/lambda)/alpha) }
  else { print("Invalid distribution") }
}