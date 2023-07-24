my_psurvreg = function(q, mean, scale = 1, distribution = "weibull", parms = NULL, lower.tail = TRUE) {
  if (!lower.tail) {
    1 - psurvreg(q = q, 
                 mean = mean, 
                 scale = scale, 
                 distribution = distribution, 
                 parms = parms)
  } else {
    psurvreg(q = q, 
             mean = mean, 
             scale = scale, 
             distribution = distribution, 
             parms = parms)
  }
}