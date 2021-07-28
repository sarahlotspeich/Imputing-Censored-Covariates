## Attempting to program the imputation method outlined in
## Atem et al (2017) using Kaplan-Meier model
# imp.model is a Kaplan-Meier model object
# surv.data is data frame containing:
#   t, the observed time min(X, C)
#   delta, the censoring indicator I(X \leq C)
csi.cox.df = function(imp.vars = c("y"), surv.data) {
  # sample size
  n = nrow(surv.data)
  
  # get imputatoin model
  imp.model <- imp.vars %>%
    # collapse vector of imputation variable names into formula RHS
    paste0(collapse = "+") %>%
    # paste LHS of formula into formula string
    paste0("survival::Surv(t, delta) ~ ", .) %>%
    # convert to formula
    as.formula() %>%
    # fit Cox model on surv.data
    survival::coxph(., data = surv.data)
  
  linear.pred <- data.matrix(surv.data[, imp.vars]) %*% matrix(data = imp.model$coeff, ncol = 1) %>%
    as.vector()
  
  # make data frame with same variable names as surv.data 
  # nit all baseline covariate values
  baseline.data = surv.data[1, ]
  baseline.data[1, ] = 0
  base.fit = survival::survfit(imp.model, baseline.data)
  
  cox.data <- data.frame(t = base.fit$time,
                         base.surv = base.fit$surv)

  surv.order <- surv.data %>%
    # exponentiate linear predictors
    mutate(exp.lin.pred = exp(linear.pred)) %>%
    # sort by observed time values
    dplyr::arrange(t) %>%
    # left join with breslow baseline cumulative hazard estimates
    left_join(cox.data, by = c("t"))
  
  # initialize imputed times as observed times
  imp.t = surv.order$t
  
  # COPIED THIS FROM SARAH: GIVE CREDIT :)
  if (surv.order[n, "delta"] == 0) {
    max.event <- max(which(surv.order[, "delta"] == 1))
    surv.max.event <- surv.order[max.event, "base.surv"]
    t.max <- surv.order[n, "t"]
    surv.order[n, "base.surv"] <- exp(t.max * log(surv.max.event))
  }
  
  t.diff = c(surv.order[-1, "t"] - surv.order[-n, "t"])
  
  for (i in which(surv.order$delta == 0)) {
    # initialize imputed value: first line of equation 5 in Atem et al (2019)
    exp.lin.pred.i = surv.order[i, "exp.lin.pred"]
    impute = (surv.order[i, "base.surv"])^(-exp.lin.pred.i)
    impute = impute/2
    
    # find terms of summation which are > 0
    positive.summands <- as.numeric(surv.order$t > surv.order[i, "t"])
    
    # baseline survival estimates raised to exponentiated linear predictor
    exp.surv.i <- (positive.summands*surv.order$base.surv)^exp.lin.pred.i
    exp.surv.sum.i <- exp.surv.i[-n] + exp.surv.i[-1]
    
    # summation on second line of equation 5
    impute = impute * sum(t.diff*exp.surv.sum.i)
    
    # add observed time according to third line of equation 5
    imp.t[i] = imp.t[i] + impute
  }
  
  # append imputed survival times to ordered survival data frame
  surv.order$x.imp = imp.t

  if(any(is.infinite(imp.t))) {print(surv.order)}
  
  return(surv.order)
}

# wrapper function for csi.cox.df which returns results from linear model fit to data
csi.cox.df.results = function(out.form = as.formula(y ~ x.imp), 
                              imp.vars = c("y"), surv.data) {
  
  # perform Cox-based conditional signle imputation
  out.model <- csi.cox.df(imp.vars = imp.vars, surv.data = surv.data) %>%
    # fit linear model to imputed data set
    lm(formula = out.form, data = .)
  
  # now, return parameter and se estimates
  return(c(coef(out.model), sqrt(diag(vcov(out.model)))))
}

## Attempting to program the imputation method outlined in
## Atem et al (2017) using Kaplan-Meier model
# surv.data is data frame containing:
#   y, the outcome for the analysis model of interest
#   t, the observed time min(X, C)
#   delta, the censoring indicator I(X \leq C)
# B is the number of bootstrap samples to be used
cmi.cox.df = function(surv.data, B = 20, p = 2,
                      out.form = as.formula(y ~ x.imp), 
                      imp.vars = c("y")) {
  n <- nrow(surv.data)  
  
  # Allocate space for simulation results
  results <- matrix(data = NA, nrow = B, ncol = 2*p)
  
  for (b in 1:B) {
    # sample with replacement to get bootstrap sample
    bs.data <- surv.data[sample(x = 1:n, size = n, replace = T), ]
    
    # Impute censored covariates of X by E[X | X > C]
    results[b, ] <- csi.cox.df.results(out.form = out.form, 
                                       imp.vars = imp.vars, 
                                       surv.data = bs.data)
  }
  
  colMeans(results) %>% return()
}