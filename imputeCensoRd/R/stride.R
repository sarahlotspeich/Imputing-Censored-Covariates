#' @import stats
pred.smooth.surv <- function(w.vector=NULL, t, data.use, covariate.value, group.value, weight)
{ 
  Xi.use = data.use[,1]
  Di.use = data.use[,2]
  Zi.use = data.use[,3]
  Wi.use = data.use[,4]
  Ui.use = data.use[,5]
  if(is.null(w.vector)) {w.vector = Wi.use}
  h = bw.nrd(Wi.use[Zi.use == covariate.value & Ui.use == group.value])
  bandwidth = h*length(Xi.use[Zi.use == covariate.value & Ui.use == group.value])^(-.10)
  K = Kern.FUN(Wi.use[Zi.use == covariate.value & Ui.use == group.value],w.vector,bandwidth)
  Xi.temp = data.use[Zi.use == covariate.value & Ui.use == group.value,1]
  Di.temp = data.use[Zi.use == covariate.value & Ui.use == group.value,2]
  Wi.temp = data.use[Zi.use == covariate.value & Ui.use == group.value,4]
  tmpind = (Xi.temp<=t)&(Di.temp==1)
  tj = Xi.temp[tmpind];
  kerni.1 = t(weight[Zi.use == covariate.value & Ui.use == group.value]*t(K))
  pihamyt0.tj.ss = sum.I(tj, "<=", Xi.temp, Vi=t(kerni.1)) ## n.tj x n.ss matrix ##
  dLamhat.tj.ss = t((kerni.1[,tmpind]))/pihamyt0.tj.ss;
  #dLamhat.tj.ss[is.na(dLamhat.tj.ss)] = 0
  ret = apply(dLamhat.tj.ss,2,sum, na.rm = TRUE)
  S.return  =exp(-ret)
  return(S.return)
}

###################################################
##SURVIVAL USING Z and W COVARIATE INFORMATION######
###################################################
#data structure: first column should be observed event time (or censoring time), second column should be indicator of whether an event was observed (Delta), third column should be discrete covariate Z, fourth column should be continuous covariate W, fifth column should indicate the U group that each person is in (for example if there are 3 "u" groups, u_1, u_2, u_3, then the fifth column would either have "1","2", or "3"); returns the data matrix with an extra column, the extra column is the survival probability at that Z and W
#newdata is an optional n by 3 matrix where the first column is the discrete covariate Z, the second column is the continuous covariate W and the third column is the U group. Predicted survival probabilities are estimated for these data;returns the data matrix with an extra column, the extra column is the survival probability at that Z and W

# in our case, t0 = 0 always

S.NPNA.zw = function(t, data, newdata = NULL, weight = NULL){
  problem <- NULL
  Xi = data[,1] # observed event/censoring time
  Di = data[,2] # censoring indicator
  Zi = data[,3] # discrete covariate
  Wi = data[,4] # continuous covariate
  Ui = data[,5] # subgroup (we will ignore or set Ui = constant for all)
  
  #if(sum(Xi > t) == 0) {print(paste("No observations past time t=",t))}
  if(is.null(weight)) {weight = rep(1,length(Xi))}
  
  zi.cat = unique(Zi)
  ui.cat = unique(Ui)
  
  ## 3/16/2017: Below will only run if newdata is NULL
  if(is.null(newdata)){
    survival.v <- rep(NA, length = dim(data)[1]) ## Changed 1/30/2017
    ##survival.v = vector(length = dim(data)[1])
    for(j in 1:length(ui.cat)) {
      ##print(j)
      for(k in 1:length(zi.cat)) {
        ##print(k)
        ui.value = ui.cat[j]
        zi.value = zi.cat[k]
        
        if(sum(Zi==zi.value & Ui ==ui.value) < 10) {
          print(paste("Warning: Very few individuals with covariate value = ",
                      zi.value, ",
                      and in U group = ",ui.value))}
        
        if(sum(Zi==zi.value & Ui == ui.value & Xi > t) < 10) {
          print(paste("Warning: Very few individuals
                      observed to survive past t=",t," with covariate value = ",
                      zi.value, ", and in U group = ",ui.value))
        }
        
        if(length( Wi[Zi == zi.value & Ui == ui.value] )<2){  ## Changed 1/30/2017
          problem <- TRUE
          print(paste("Warning-Error: Less than 2 individuals with W covariate for z-covariate value = ",
                      zi.value, ", and in U group = ",ui.value,
                      ". Setting estimate=0.  W is",
                      Wi[Zi == zi.value & Ui == ui.value ], sep=""))
          P.return <- 0
        } else {
          P.return <- pred.smooth.surv(w.vector = Wi[Zi == zi.value & Ui == ui.value], t=t,
                                       data.use = data, covariate.value = 	zi.value,
                                       group.value = ui.value, weight = weight)
        }
        survival.v[Zi == zi.value & Ui == ui.value] <- P.return
      }
    }
    data = cbind(data, survival.v)
  }
  
  if(!is.null(newdata)) {
    Zi.new = newdata[,1]
    Wi.new = newdata[,2]
    Ui.new = newdata[,3]
    
    survival.v.new = rep(NA,length = dim(newdata)[1]) ## Changed 1/30/2017
    ##survival.v.new = vector(length = dim(newdata)[1])
    for(j in 1:length(ui.cat)) {
      for(k in 1:length(zi.cat)) {
        ui.value = ui.cat[j]
        zi.value = zi.cat[k]
        
        ##3/16/2017: Added warning message
        #if(sum(Zi==zi.value & Ui ==ui.value) < 10) {
        #  print(paste("Warning: Very few individuals with covariate value = ",
        #  			zi.value, ",
        #    and in U group = ",ui.value))}
        
        ##3/16/2017: Added warning message
        #if(sum(Zi==zi.value & Ui == ui.value & Xi > t) < 10) {
        #  print(paste("Warning: Very few individuals
        #                  observed to survive past t=",t," with covariate value = ",
        #                  zi.value, ", and in U group = ",ui.value))
        #}
        
        
        if(length( Wi[Zi == zi.value & Ui == ui.value] )<2){ ## Changed 1/30/2017
          problem <- TRUE
          print(paste("Warning-Error: Less than 2 individuals with W covariate
                      for z-covariate value = ",
                      zi.value, ", and in U group = ",ui.value,
                      ". Setting estimate=0. W is",
                      Wi[Zi == zi.value & Ui == ui.value],
                      sep=""))
          P.return <- 0
        } else {
          P.return = pred.smooth.surv(w.vector = Wi.new[Zi.new == zi.value &
                                                          Ui.new == ui.value],
                                      t=t,data.use = data, covariate.value =  zi.value, group.value = ui.value,
                                      weight = weight)
        }
        survival.v.new[Zi.new == zi.value & Ui.new == ui.value] = P.return
      }
    }
    newdata.matrix = cbind(newdata, survival.v.new)
  }
  if(is.null(newdata)) {return(list(data.out=data,problem=problem))}
  if(!is.null(newdata)) {return(list(data.out=newdata.matrix,problem=problem))}
}
