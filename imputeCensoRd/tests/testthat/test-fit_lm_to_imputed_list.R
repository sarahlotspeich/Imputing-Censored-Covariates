# 4 tests in total
set.seed(95)
library(usethis)

# tests that bad input produces simple errors
# 3 tests in total
test_that("simple errors for bad input", {
  # generate sample data and KM fit for test
  sample.data <- generate_data(n = 10, n.sims = 1, beta0 = 0, betaX = 1)
  # perform conditional mean imputation with bootstrap
  imp.sample.data <- condl_mean_impute_bootstrap(obs = "t", event = "event",
                                                 data = sample.data, M = 20)
  
  # list contains elements that is neither a dataframe nor a matrix
  bad.sample.data <- imp.sample.data
  bad.sample.data[[20]] <- "character"
  expect_error(fit_lm_to_imputed_list(imputed_list = bad.sample.data, formula = as.formula(y ~ x)), "data.frame or")
  
  # one dataframe in list has a misnamed column
  bad.sample.data <- imp.sample.data
  colnames(bad.sample.data[[20]])[4] <- "Y"
  expect_error(fit_lm_to_imputed_list(imputed_list = bad.sample.data, formula = as.formula(y ~ x)), "needed for formula")
  
  # formula contains variable not found in any list
  expect_error(fit_lm_to_imputed_list(imputed_list = imp.sample.data, formula = as.formula(y ~ s)), "needed for formula")
})

# test that the returned list has expected properties
# 1 test in total
test_that("test for proper output", {
  # generate sample data and KM fit for test
  sample.data <- generate_data(n = 10, n.sims = 1, beta0 = 0, betaX = 1)
  sample.fit <- with(sample.data, survival::survfit(formula = survival::Surv(x, event) ~ 1))
  # perform conditional mean imputation
  imp.sample.data <- condl_mean_impute_bootstrap(obs = "t", event = "event",
                                                 data = sample.data, M = 20)
  est.summary <- fit_lm_to_imputed_list(imputed_list = imp.sample.data, formula = as.formula(y ~ imp))
  
  # check that the variance estimates obey simple rule
  expect_true(all(est.summary$Pooled_Var >= 0))
})