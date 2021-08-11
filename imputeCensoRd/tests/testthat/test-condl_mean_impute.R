# 13 tests in total
set.seed(95)
library(usethis)

# tests that bad input produces simple errors
# 12 tests in total
test_that("simple errors for bad input", {
  # generate sample data and KM fit
  sample.data <- generate_data(n = 10, n.sims = 1, beta0 = 0, betaX = 1)
  sample.fit <- with(sample.data, survival::survfit(formula = survival::Surv(x, event) ~ 1))
  
  # no arguments given
  expect_error()

  # time, event, or addl_covar not characters
  expect_error(condl_mean_impute(fit = sample.fit, obs = 2, event = "event", data = sample.data), "obs must be a character")
  expect_error(condl_mean_impute(fit = sample.fit, obs = "t", event = 2, data = sample.data), "event must be a character")
  expect_error(condl_mean_impute(fit = sample.fit, obs = "t", event = "event", addl_covar = 2, data = sample.data), "addl_covar must be a character")
  
  # approx_beyond not an acceptable option
  expect_error(condl_mean_impute(fit = sample.fit, obs = "t", event = "event", data = sample.data, approx_beyond = "none"), "expo, zero, or carryforward")
  
  # data not a data frame or a matrix
  expect_error(condl_mean_impute(fit = sample.fit, obs = "t", event = "event", data = 2), "data must be a data frame")
  
  # data does not contain specified columns
  expect_error(condl_mean_impute(fit = sample.fit, obs = "obs", event = "event", data = sample.data), "column with name obs")
  expect_error(condl_mean_impute(fit = sample.fit, obs = "t", event = "d", data = sample.data), "column with name d")
  expect_error(condl_mean_impute(fit = sample.fit, obs = "t", event = "event", addl_covar = "a", data = sample.data), "columns with names: a")
  expect_error(condl_mean_impute(fit = sample.fit, obs = "t", event = "event", addl_covar = c("a", "z"), data = sample.data), "columns with names: a, z")
  
  # columns of data have values outside of expected range
  # negative t
  bad.sample.data = sample.data
  bad.sample.data[1, "t"] = -1
  expect_warning(condl_mean_impute(fit = sample.fit, obs = "t", event = "event", data = bad.sample.data), "must be positive")
  # event not 0 or 1
  bad.sample.data = sample.data
  bad.sample.data[1, "event"] = 2
  expect_warning(condl_mean_impute(fit = sample.fit, obs = "t", event = "event", data = bad.sample.data), "must be either 0 or 1")
})

# test that elements of the returned data have expected properties
# 1 test in total
test_that("test for proper output", {
  # generate sample data and KM fit for test
  sample.data <- generate_data(n = 10, n.sims = 1, beta0 = 0, betaX = 1)
  sample.fit <- with(sample.data, survival::survfit(formula = survival::Surv(x, event) ~ 1))
  # perform conditional mean imputation
  imp.sample.data <- condl_mean_impute(fit = sample.fit, obs = "t", event = "event", data = sample.data)
  
  # expect that imputed values are greater than or equal to observed times
  expect_true(all(imp.sample.data$t <= imp.sample.data$imp))
})