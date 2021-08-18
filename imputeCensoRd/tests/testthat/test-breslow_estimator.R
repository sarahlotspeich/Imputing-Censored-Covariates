set.seed(95)
library(usethis)

# tests that bad input produces simple errors
# 11 tests in total
test_that("simple errors for bad input", {
  # no arguments given
  expect_error(breslow_estimator())
  
  # generate sample data and fit imputation model
  sample.data <- generate_data(n = 10, n.sims = 1, beta0 = 0, betaX = 1)
  sample.fit <- survival::coxph(formula = survival::Surv(t, event) ~ y, data = sample.data)
  sample.data$hr <- exp(sample.fit$linear.predictors)

  # time, event, or hr not characters
  expect_error(breslow_estimator(time = 2, event = "event", hr = "hr", data = sample.data), "time must be a character")
  expect_error(breslow_estimator(time = "t", event = 2, hr = "hr", data = sample.data), "event must be a character")
  expect_error(breslow_estimator(time = "t", event = "event", hr = 2, data = sample.data), "hr must be a character")
  
  # data not a data frame
  expect_error(breslow_estimator(time = "t", event = "event", hr = "hr", data = 2), "data must be a data frame")
  
  # data does not contain specified columns
  expect_error(breslow_estimator(time = "time", event = "event", hr = "hr", data = sample.data), "column with name time")
  expect_error(breslow_estimator(time = "t", event = "d", hr = "hr", data = sample.data), "column with name d")
  expect_error(breslow_estimator(time = "t", event = "event", hr = "HazardRatio", data = sample.data), "column with name HazardRatio")
  
  # columns of data have values outside of expected range
  # negative t
  bad.sample.data = sample.data
  bad.sample.data[1, "t"] = -1
  expect_warning(breslow_estimator(time = "t", event = "event", hr = "hr", data = bad.sample.data), "must be positive")
  # event not 0 or 1
  bad.sample.data = sample.data
  bad.sample.data[1, "event"] = 2
  expect_warning(breslow_estimator(time = "t", event = "event", hr = "hr", data = bad.sample.data), "either 0 or 1")
  # negative hr
  bad.sample.data = sample.data
  bad.sample.data[1, "hr"] = -1
  expect_warning(breslow_estimator(time = "t", event = "event", hr = "hr", data = bad.sample.data), "between 0 and 1")
})

# test that elements of the returned have expected properties
# 2 tests in total
test_that("test for proper output", {
  # generate sample data and fit imputation model
  sample.data <- generate_data(n = 10, n.sims = 1, beta0 = 0, betaX = 1)
  sample.fit <- survival::coxph(formula = survival::Surv(t, event) ~ y, data = sample.data)
  sample.data$hr <- exp(sample.fit$linear.predictors)

  # calculate breslow estimator
  sample.breslow <- breslow_estimator(time = "t", event = "event", hr = "hr", data = sample.data)
  # times element from breslow_estimator should match the sorted, unique event times in sample.data
  expect_true(all(sample.breslow$times == sort(unique(filter(sample.data, event == 1)$t))))
  # baseline survival estimates should all be inclusively between 0 and 1
  expect_true(all(sample.breslow$basesurv <= 1 & sample.breslow$basesurv >= 0))
})
