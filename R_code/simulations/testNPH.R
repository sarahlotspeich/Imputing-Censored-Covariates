rm(list = ls())

# Sarah's directory
# source("~sarahlotspeich/Dropbox/UNC/Sarah-Lotspeich/code/mi/impute_multivar.R")

library(tidyverse)

# Kyle's directory
setwd("~kylegrosser/Documents/GitHub/Imputing-Censored-Covariates")
source("R_code/source/Generate_Data.R")
source("R_code/source/Cox_Imputation_DF.R")

sim_seed <- 918
set.seed(sim_seed)

PH <- FALSE
N <- 500
beta0 <- 1
beta1 <- 1
beta2 <- 0.25
pZ <- 0.5
n.sims <- 50

res <- data.frame()
for (scaleC in c(0.35)) {
    new_res <- data.frame(sim = 1:n.sims, scaleC, N, PH, 
                          b0 = NA, b1 = NA, b2 = NA, 
                          b0_se = NA, b1_se = NA, b2_se = NA)
    for (s in 1:n.sims) {
        # Simulate data
        z <- rbinom(n = N, size = 1, prob = pZ)
        if (PH) {
          # Proportional hazard for X
          x <- rweibull(n = N, shape = 0.75, scale = 0.25)
        } else {
          # Non-proportional hazard for X
          x <- rweibull(n = N, shape = (0.6 + 0.9 * z), scale = 0.25)
        }
        e <- rnorm(n = N, mean = 0, sd = 1)
        y <- beta0 + beta1 * x + beta2 * z + e
        c <- rweibull(n = N, shape = 1, scale = scaleC) # for light or scale = 0.35 for heavy 
        delta <- as.numeric(x <= c)
        t <- pmin(x, c)
        x[delta == 0] <- NA
        sim_dat <- data.frame(y, x, t, z, delta)
        # CMI - Atem et al. (2017)
        ## Using survfit() to estimate baseline survival
        cmi_sf <- cmi.cox.df(surv.data = sim_dat, B = 20, p = 3, 
                             out.form = as.formula("y ~ x.imp + z"),
                             imp.vars = c("y", "z"))
        new_res[s, c("b0", "b1", "b2")] <- cmi_sf[1:3]
        new_res[s, c("b0_se", "b1_se", "b2_se")] <- cmi_sf[4:6]
    }
    res <- rbind(res, new_res)
}

res[5:10] %>% colMeans()
