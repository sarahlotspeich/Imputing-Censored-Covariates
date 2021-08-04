test_that("simple errors for bad input", {
  # no arguments given
  expect_error(generate_data())
  
  # n or n.sims not positive integers
  expect_error(generate_data(n = -1, beta0 = 0, betaX = 1))
  expect_error(generate_data(n = 1.5, beta0 = 0, betaX = 1))
  
  # invalid parameters for data simulation
  expect_error(generate_data(n = 10, n.sims = -1, beta0 = 0, betaX = 1))
  expect_error(generate_data(n = 10, beta0 = 0, betaX = 1, xshape = 0))
  expect_error(generate_data(n = 10, beta0 = 0, betaX = 1, c.lower = 5, c.upper = 1))
  
  # number of observed covariates does not match number of corresponding parameters
  betaZ = c(1.5)
  Z = matrix(1, 10, 2)
  expect_error(generate_data(n = 10, beta0 = 0, betaX = 1, betaZ = betaZ, z = z))
})
