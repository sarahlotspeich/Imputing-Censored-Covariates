# Correcting conditional mean imputation for censored covariates and improving usability
## Lotspeich, Grosser & Garcia 
### Statistical methods to impute censored covariates. 

The complete R package `imputeCensoRd` and code for the simulations included in the note.

# Repo Organization 

- bash_scripts: Bash scripts used to run R simulations on cluster
- R_code: contains R code for
	- markdowns: RMD files used to show results from simulations
	- simulations: R scripts to run simulations of censored covariate imputation
	- source: functions used in simulations
- Simulation_results: contain RDS files created by simulations

# Example

![](Sim-Setup.png)

```{r}
sim_seed <- 114

# Set parameters 
N <- 1000
lambda <- -2
beta0 <- 1
beta1 <- 1
beta2 <- 0.25

# Simulate data
z <- rbinom(n = N, size = 1, prob = 0.25)
x <- sim_cox(n = N, logHR = lambda, covariate = matrix(z, ncol = 1), dist = "Exponential", lambda = 5)
e <- rnorm(n = N, mean = 0, sd = 1)
y <- beta0 + beta1 * x + beta2 * z + e
c <- rexp(n = N, rate = 4)
delta <- as.numeric(x <= c)
t <- pmin(x, c)
x[delta == 0] <- NA
sim_dat <- data.frame(y, x, t, z, delta)
```

# References

Bender, R., Augustin, T., and Blettner, M. (2005). Generating survival times to simulate Cox proportional hazards models. *Statistics in Medicine*, 24:1713–1723.
