source("atem_impute_multivar.R")

sim_cox <- function(n, logHR, covariate, dist = "Weibull", lambda = 0.25, nu = 0.75, alpha = NULL) {
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

args <- commandArgs(TRUE)
sim_seed <- 114 + as.integer(args)

pZ.to.test <- c(0.25)
logHR.to.test <- seq(-2, 2)
N.to.test <- c(100, 500, 1000)

scaleC <- 2
beta0 <- 1
beta1 <- 1
beta2 <- 0.25

n.sims <- 50

for (N in N.to.test) {
  for (logHR in logHR.to.test) {
    res <- data.frame()
    set.seed(sim_seed)
    for (pZ in pZ.to.test) {
        for (s in 1:n.sims) {
        # Simulate data
        z <- rbinom(n = N, size = 1, prob = pZ)
        x <- sim_cox(n = N, logHR = logHR, covariate = matrix(z, ncol = 1), dist = "Exponential", lambda = 5)
        #sim_cox(n = N, logHR = logHR, covariate = matrix(z, ncol = 1), dist = "Weibull", lambda = 1, nu = 1)
        e <- rnorm(n = N, mean = 0, sd = 1)
        y <- beta0 + beta1 * x + beta2 * z + e
        c <- rexp(n = N, rate = 4)
        #c <- rweibull(n = N, shape = 1, scale = scaleC)
        delta <- as.numeric(x <= c)
        t <- pmin(x, c)
        x[delta == 0] <- NA
        sim_dat <- data.frame(y, x, t, z, delta)
        
        for (form in c("JRSS-C", "SMMR", "BIOMJ", "CORRECT")) {
          fit <- atem_cmi_multivar(data = sim_dat, outcome = "y", time = "t", event = "delta", covar = "z", 
                                    approx_last = "expo", est_S0 = "breslow", atem_formula = form, M = 20, indicator_incl = FALSE)#, 
                                    #save_imputed_as = paste0(ifelse(form == "JRSS-C", "JRSSC", form), "_pZ", pZ * 100, "_N", N, "_logHR", ifelse(logHR < 0, "neg", ""), abs(logHR), "_sim", s, "_seed", sim_seed, ".csv"))
          res <- rbind(res, 
                       data.frame(sim = s, N, pZ, logHR, cens = mean(1 - delta), atem_formula = form,
                                  b0 = fit$coeff[1], b1 = fit$coeff[2], b2 = fit$coeff[3], 
                                  b0_se = fit$se[1], b1_se = fit$se[2], b2_se = fit$se[3]))
        }
        write.csv(res, paste0("EXCL_vary_logHR", ifelse(logHR < 0, "neg", ""), abs(logHR), "_n", N, "_seed", sim_seed, ".csv"), row.names = F)
      }
    }
  } 
}
