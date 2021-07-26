### Clear workspace
rm(list = ls())

print("Scenario 1: Light censoring")

### Dependencies
library(tidyverse)
library(data.table)
library(survival)

# setwd("Documents/GitHub/re-impute/Misspecification_Studies/Experiments_Set1")
### R scripts
source("R_scripts/Generate_Data.R")
source("R_scripts/Cox_Imputation_DF.R")
source("R_scripts/KM_Imputation.R")

### Simulation parameters
# Sample size
n.to.test <- c(500)
# Imputation models
models.to.test <- list(c("y", "a"), c("y"), c("a"), "KM")
# Number of simulations
n.sims <- 100
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

# base.data = sim.data[1, ]
# base.data[1, ] = 0
# mod = coxph(Surv(t, delta) ~ y + a, sim.data)
# base.data
# base.fit = survival::survfit(mod, base.data)
# base.data = data.frame(base.surv = base.fit$surv,
#                        t = base.fit$time)
# with.base.surv = sim.data %>%
#   left_join(base.data, by = "t")


sim.result = list(Detail = paste("n.sims =", n.sims, "censoring rate =", 1 - mean(sim.data$delta)),
                  Table = results.tab)

saveRDS(sim.result, "result_tables/Scenario_1_Results.RDS")
readRDS("result_tables/Scenario_1_Results.RDS")