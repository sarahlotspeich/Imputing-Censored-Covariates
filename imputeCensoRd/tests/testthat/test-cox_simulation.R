# tests that bad input produces simple errors
# 4 tests in total
test_that("simple errors for bad input", {
  # no arguments given
  expect_error(cox_simulation())
  
  n = 10
  logHR = 1
  A = matrix(rep(c(0, 1), 5), n)
  
  # n not a positive integer
  expect_error(cox_simulation(n = -1, logHR = logHR, A = A), "positive integer")
  expect_error(cox_simulation(n = 1.5, logHR = logHR, A = A), "positive integer")
  
  # invalid distribution character
  expect_error(cox_simulation(n = n, logHR = logHR, A = A, dist = "nope"), "Exponential, Weibull, or Gompertz")
  
  # number of observed covariates does not match number of corresponding parameters
  lorHR = c(1, 1.5)
  expect_error(cox_simulation(n = n, logHR = lorHR, A = A), "length")
})