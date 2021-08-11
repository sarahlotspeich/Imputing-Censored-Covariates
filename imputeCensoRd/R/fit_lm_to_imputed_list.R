fit_lm_to_imputed_list = function(imputed.list, formula) {
  # determine number of bootstrap samples
  M = length(imputed.list)
  
  # get list of length M, with each entry a vector of coefficient and SE estimates
  estimates.list <- lapply(X = imp.sample.data, 
                           FUN = function(x) {
                             model <- lm(formula = formula, data = x)
                             estimates <- c(model$coefficients, sqrt(diag(vcov(model))))
                           })

  # collapse list of parameter estimates into M x p matrix
  estimates.matrix <- matrix(0, nrow = M, ncol = length(param.estimates[[1]]))
  for (i in 1:M) estimates.matrix[i, ] <- estimates.list[[i]]
  
  # return column means
  colMeans(estimates.matrix)
}

# # generate sample data and KM fit for test
# sample.data <- generate_data(n = 500, n.sims = 1, beta0 = 0, betaX = 1)
# sample.fit <- with(sample.data, survival::survfit(formula = survival::Surv(x, event) ~ 1))
# # perform conditional mean imputation
# imp.sample.data <- condl_mean_impute_bootstrap(fit = sample.fit, obs = "t", event = "event", 
#                                                data = sample.data, M = 20)
# fit_to_imputed_data(imputed.list = imp.sample.data, formula = as.formula(y ~ imp))