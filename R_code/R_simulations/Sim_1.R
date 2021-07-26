### Clear workspace
rm(list = ls())

print("Simulation 1: One Cox covariate with one auxiliary variable; Light censoring")

### Dependencies
library(tidyverse)
library(survival)

# setwd("../../")
### R scripts
getwd()
source("R_code/R_source/Generate_Data.R")
source("R_code/R_source/Cox_Imputation_DF.R")
source("R_code/R_source/KM_Imputation.R")

### Simulation parameters
# Sample size
n.to.test <- c(100, 500, 1000)
# Imputation models
models.to.test <- list(c("y", "a"), c("y"), c("a"), "KM")
# Number of simulations
n.sims <- 1000
# Outcome model parameters
beta0 <- 1; beta1 <- 0.5
# Regression, shape, and scale parameters in TRUE model for generation of X
gamma1 <- 1
xscale <- 0.5
xshape <- 1.5
# Censoring mechanism (C) parameter for ~20% censoring
c.upper <- 5.5
# Number of bootstrap samples
B <- 20

### Simulation with Auxiliary Variables
set.seed(95)

all.results = data.frame()
for (n in n.to.test) {
  # Total number of samples
  N <- n*n.sims

  # Generate auxiliary random variables for Cox simuilation
  A = rbinom(n = N, size = 1, prob = 0.25) %>%
    matrix(N, 1)

  # Generate covariate x for outcome model based on Cox simulation
  x = cox.simulation(n = N, param = gamma1, covariate = A, 
                     dist = "Weibull", lambda = xscale, nu = xshape)

  # Generate data ONCE to save time in simulation
  sim.data <- generate.data(n = n, n.sims = n.sims, 
                            beta0 = beta0, betaX = beta1, 
                            x = x, c.upper = c.upper)
  
  sim.data$a <- A
  sim.data <- sim.data %>%
    data.frame()

  n.results = data.frame()
  for (model in models.to.test) {

    # Allocate space for simulation results
    results <- data.frame(rep = 1:n.sims, 
                          beta0.est = NA, beta1.est = NA, 
                          beta0.se = NA, beta1.se = NA,
                          imputation.model = paste0(model, collapse = "+"),
                          n = n)

    # Run simulation for sample size n
    for (sim in 1:n.sims) {
      # Grab simulation data by simulation id (sim.id)
      model.data <- sim.data %>%
        dplyr::filter(sim.id == sim)
  
      # Store model parameter and SE estimates
      if (any(model == "KM")) {
        results[sim, 2:5] = cmi.km(surv.data = model.data, B = B)
      }
      else {
        results[sim, 2:5] = cmi.cox.df(surv.data = model.data, B = B,
                                       imp.vars = model)
      }
    }
    
    n.results = rbind(n.results, results)
    print(paste("Sim. with model", paste0(model, collapse = "+"), "and n =", n, "complete"))
  }
  all.results = rbind(all.results, n.results)
}

results.tab <- all.results %>% 
  group_by(n, imputation.model) %>%
  summarize(beta0.estimate = mean(beta0.est),
            beta0.se = mean(beta0.se),
            beta1.estimate = mean(beta1.est),
            beta1.se = mean(beta1.se))

sim.result = list(Detail = paste("n.sims =", n.sims, "censoring rate =", 1 - mean(sim.data$delta)),
                  Table = results.tab)

saveRDS(sim.result, "Simulation_Results/Sim_1_Results.RDS")
readRDS("Simulation_Results/Sim_1_Results.RDS")