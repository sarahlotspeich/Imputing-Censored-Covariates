# 13 tests in total
set.seed(95)
library(usethis)

# tests that bad input produces simple errors
# 10 tests in total
test_that("simple errors for bad input", {
  # generate sample data and KM fit
  sample.data <- generate_data(n = 10, n.sims = 1, beta0 = 0, betaX = 1)

  # no arguments given
  expect_error(condl_mean_impute_bootstrap())
  
  # time, event, or addl_covar not characters
  expect_error(condl_mean_impute_bootstrap(obs = 2, event = "event", data = sample.data, M = 20), "obs must be a character")
  expect_error(condl_mean_impute_bootstrap(obs = "t", event = 2, data = sample.data, M = 20), "event must be a character")
  expect_error(condl_mean_impute_bootstrap(obs = "t", event = "event", addl_covar = 2, data = sample.data, M = 20), "addl_covar must be a character")
  
  # approx_beyond not an acceptable option
  expect_error(condl_mean_impute_bootstrap(obs = "t", event = "event", data = sample.data, approx_beyond = "none", M = 20), "expo, zero, or carryforward")
  
  # data not a data frame or a matrix
  expect_error(condl_mean_impute_bootstrap(obs = "t", event = "event", data = 2, M = 20), "data must be a data frame")
  
  # data does not contain specified columns
  expect_error(condl_mean_impute_bootstrap(obs = "obs", event = "event", data = sample.data, M = 20), "column with name obs")
  expect_error(condl_mean_impute_bootstrap(obs = "t", event = "d", data = sample.data, M = 20), "column with name d")
  expect_error(condl_mean_impute_bootstrap(obs = "t", event = "event", addl_covar = "a", data = sample.data, M = 20), "columns with names: a")
  expect_error(condl_mean_impute_bootstrap(obs = "t", event = "event", addl_covar = c("a", "z"), data = sample.data, M = 20), "columns with names: a, z")
})

# test that the returned list has expected properties
# 3 test in total
test_that("test for proper output", {
  # generate sample data and KM fit for test
  sample.data <- generate_data(n = 10, n.sims = 1, beta0 = 0, betaX = 1)
  sample.fit <- with(sample.data, survival::survfit(formula = survival::Surv(x, event) ~ 1))
  # perform conditional mean imputation
  imp.sample.data <- condl_mean_impute_bootstrap(obs = "t", event = "event", 
                                                 data = sample.data, M = 20)

  # expect 20 data frames in returned list
  expect_true(length(imp.sample.data) == 20)
  expect_true(all(unlist(lapply(X = imp.sample.data, is.data.frame))))
  # expect that imputed values are greater than or equal to observed times
  expect_true(all(unlist(lapply(X = imp.sample.data, FUN = function(x) all(x[, "t"] <= x[, "imp"])))))
})