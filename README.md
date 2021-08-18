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

## Setup 

Samples of `n=1000` subjects were created, beginning with fully-observed binary covariate `Z` which was generated from a Bernoulli distribution with P(Z = 1) = 0.25. The right censoring variable, `C`, was independently generated from its own exponential distribution with rate=4. Covariate `X` was generated to have an exponential baseline survival function with rate=5 following the procedure of Bender et al. (2005). To do so, we generated U from a uniform distribution with min=0 and max=1 and constructed the baseline cumulative hazard as $H_0(X) = -\log(U)\exp(-\lambda Z)$, from which we have $X = H_0(X) / \lambda$. The $\lambda$ used to generate `X` is the log hazard ratio, and values of $\lambda = -2, -1, 0, 1, 2$ were considered, leading to an average of 55%, 51%, 45%, 39%, and 0.36% censoring, respectively. Finally, observed values were constructed as T = min(X, C). The outcome, $Y$, is calculated as `Y = 1 + X + 0.25Z + e`, with random errors `e` generated independently from a standard normal distribution. Each setting was replicated 1000 times. 

# References

Bender, R., Augustin, T., and Blettner, M. (2005). Generating survival times to simulate cox proportional hazardsmodels. *Statistics in Medicine*, 24:1713–1723.
