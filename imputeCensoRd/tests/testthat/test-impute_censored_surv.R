set.seed(95)
library(usethis)

# tests that bad input produces simple errors
# 11 tests in total
test_that("simple errors for bad input", {
  # no arguments given
  expect_error(impute_censored_surv())
  
  # generate sample data and fit imputation model
  sample.data <- generate_data(n = 10, n.sims = 1, beta0 = 0, betaX = 1)
  sample.fit <- survival::coxph(formula = survival::Surv(t, event) ~ y, data = sample.data)
  sample.data$hr <- exp(sample.fit$linear.predictors)
  
  # calculate breslow estimator
  sample.breslow <- breslow_estimator(time = "t", event = "event", hr = "hr", data = sample.data)
  surv.data <- with(sample.breslow, data.frame(t = times, base.surv = basesurv))
  sample.data <- sample.data %>%
    dplyr::left_join(surv.data, by = "t")
  at_time <- subset(sample.data, event == 1)[1, "t"]

  # time, event, or surv not characters
  expect_error(impute_censored_surv(at_time = at_time, time = 2, event = "event", surv = "base.surv", data = sample.data), "time must be a character")
  expect_error(impute_censored_surv(at_time = at_time, time = "t", event = 2, surv = "base.surv", data = sample.data), "event must be a character")
  expect_error(impute_censored_surv(at_time = at_time, time = "t", event = "event", surv = 2, data = sample.data), "surv must be a character")
  
  # data not a data frame
  expect_error(impute_censored_surv(at_time = at_time, time = "t", event = "event", surv = "base.surv", data = 2), "data must be a data frame")
  
  # data does not contain specified columns
  expect_error(impute_censored_surv(at_time = at_time, time = "time", event = "event", surv = "base.surv", data = sample.data), "column with name time")
  expect_error(impute_censored_surv(at_time = at_time, time = "t", event = "d", surv = "base.surv", data = sample.data), "column with name d")
  expect_error(impute_censored_surv(at_time = at_time, time = "t", event = "event", surv = "surv", data = sample.data), "column with name surv")
  
  # columns of data have values outside of expected range
  # negative t
  bad.sample.data = sample.data
  bad.sample.data[1, "t"] = -1
  bad.sample.data
  expect_warning(impute_censored_surv(at_time = at_time, time = "t", event = "event", surv = "base.surv", data = bad.sample.data), "must be positive")
  # event not 0 or 1
  bad.sample.data = sample.data
  bad.sample.data[1, "event"] = 2
  expect_warning(impute_censored_surv(at_time = at_time, time = "t", event = "event", surv = "base.surv", data = bad.sample.data), "either 0 or 1")
  # negative base.surv
  bad.sample.data = sample.data
  bad.sample.data[1, "base.surv"] = -1
  expect_warning(impute_censored_surv(at_time = at_time, time = "t", event = "event", surv = "base.surv", data = bad.sample.data), "between 0 and 1")
})
