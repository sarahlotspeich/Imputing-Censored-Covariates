## Attempting to program the imputation method outlined in
## Atem et al (2017) using Kaplan-Meier model
# surv.data is data frame containing:
#   t, the observed time min(X, C)
#   delta, the censoring indicator I(X \leq C)
csi.km = function(surv.data) {
  # sample size
  n = nrow(surv.data)
  
  # calculate KM estimate
  imp.model = survfit(survival::Surv(t, delta) ~ 1, data = surv.data)
  
  # get KM estimate data
  km.data <- data.frame(t = imp.model$time, 
                        delta = imp.model$n.event, 
                        surv = imp.model$surv)
  
  surv.order <- surv.data %>%
    # sort by increasing values of t
    arrange(t)
  
  # capture survival function estimate
  surv.order = suppressMessages(surv.order <- left_join(surv.order, km.data, join.by = c("t", "delta")))
  
  # initialize imputed times as observed times
  imp.t = surv.order$t
  
  # COPIED THIS FROM SARAH: GIVE CREDIT :)
  if (surv.order[n, "delta"] == 0) {
    max.event <- max(which(surv.order[, "delta"] == 1))
    surv.max.event <- surv.order[max.event, "surv"]
    t.max <- surv.order[nrow(surv.order), "t"]
    surv.order[n, "surv"] <- exp(t.max * log(surv.max.event))
  }
  
  surv.sum = surv.order[1:(n - 1), "surv"] + surv.order[2:n, "surv"]
  t.diff = surv.order[2:n, "t"] - surv.order[1:(n - 1), "t"]
  surv.sum.times.t.diff = c(surv.sum*t.diff, 0)
  
  for (i in which(surv.order$delta == 0)) {
    # initialize imputed value: first line of equation 5 in Atem et al (2019)
    impute = (2*surv.order[i, "surv"])^(-1)
    
    # summation on second line of equation 5
    impute = impute * sum(surv.sum.times.t.diff[which(surv.order$t > surv.order[i, "t"])])
    
    # add observed time according to third line of equation 5
    imp.t[i] = imp.t[i] + impute
  }
  
  # append imputed survival times to ordered survival data frame
  surv.order$x.imp = imp.t
  
  return(surv.order)
}

# wrapper function for csi.km which returns results from linear model fit to data
csi.km.results = function(out.form = as.formula(y ~ x.imp),
                          surv.data) {
  # perform KM-based conditional single imputation
  out.model <- csi.km(surv.data = surv.data) %>%
    # fit linear model to imputed data set
    lm(formula = out.form, data = .)
  
  # now, return parameter and se estimates
  return(c(coef(out.model), sqrt(diag(vcov(out.model)))))
}

## Attempting to program the multiple imputation method outlined in
## Atem et al (2017) using Kaplan-Meier model
# surv.data is data frame containing:
#   y, the outcome for the analysis model of interest
#   t, the observed time min(X, C)
#   delta, the censoring indicator I(X \leq C)
# B is the number of bootstrap samples to be used
cmi.km = function(surv.data, B = 20, p = 2,
                  out.form = as.formula(y ~ x.imp)) {
  n <- nrow(surv.data)  
  
  # Allocate space for simulation results
  results <- matrix(data = NA, nrow = B, ncol = 2*p)
  
  for (b in 1:B) {
    # sample with replacement to get bootstrap sample
    bs.data <- surv.data[sample(x = 1:n, size = n, replace = T), ]
    
    # Impute censored covariates of X by E[X | X > C]
    results[b, ] <- csi.km.results(surv.data = bs.data, out.form = out.form)
  }
  
  colMeans(results) %>% return()
}