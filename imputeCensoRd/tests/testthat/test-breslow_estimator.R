# tests that bad input produces simple errors
# 6 tests in total
test_that("simple errors for bad input", {
  # no arguments given
  expect_error(breslow_estimator())
  
  # # n not a positive integer
  # expect_error(generate_data(n = -1, beta0 = 0, betaX = 1), "positive integer")
  # expect_error(generate_data(n = 1.5, beta0 = 0, betaX = 1), "positive integer")
  # 
  # # invalid parameters for data simulation
  # expect_error(generate_data(n = 10, beta0 = 0, betaX = 1, xshape = 0), "must be positive")
  # expect_error(generate_data(n = 10, beta0 = 0, betaX = 1, c.lower = 5, c.upper = 1), "must be less than")
  # 
  # # number of observed covariates does not match number of corresponding parameters
  # betaZ = c(1.5)
  # z = matrix(1, 10, 2)
  # expect_error(generate_data(n = 10, beta0 = 0, betaX = 1, betaZ = betaZ, z = z), "must equal number of columns")
})

# # test that data frame returned is proper
# # 5 tests in total
# test_that("test for proper output", {
#   # generate data without observed covariates
#   sample.result.noZ = generate_data(n = 10, beta0 = 0, betaX = 1)
#   betaZ = c(1.5, 0.5)
#   z = matrix(1, 10, 2)
#   
#   # generate data with observed covariates
#   sample.result.withZ = generate_data(n = 10, beta0 = 0, betaX = 1, betaZ = betaZ, z = z)
#   
#   # dimensions are correct
#   expect_true(nrow(sample.result.noZ) == 10)
#   expect_true(ncol(sample.result.noZ) == 6)
#   expect_true(ncol(sample.result.withZ) == 8)
#   
#   # t is indeed less than or equal to both x and c
#   expect_true(all(sample.result.noZ$t <= sample.result.noZ$x))
#   expect_true(all(sample.result.noZ$t <= sample.result.noZ$c))
# })