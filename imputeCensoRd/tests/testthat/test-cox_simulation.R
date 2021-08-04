# tests that bad input produces simple errors
# 4 tests in total
test_that("simple errors for bad input", {
  # no arguments given
  expect_error(cox_simulation())
  
  logHR = 1
  A = rep(c(0, 1), 5)
  
  # n not a positive integer
  expect_error(cox_simulation(n = -1, logHR = logHR, A = A))
  expect_error(cox_simulation(n = 1.5, logHR = logHR, A = A))
  
  # invalid distribution character
  expect_error(cox_simulation(n = 10, logHR = logHR, A = A, dist = "nope"))
})