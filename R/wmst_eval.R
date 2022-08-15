################################
# wmst_eval (two-arm) -- hidden
################################

wmst_eval<-function(survC, survT, eperiod, end_study, tau0, tau1, n){
  if (tau1<=end_study - eperiod) return(1)
  P_x_followup<-function(x){ # density for exactly x amount of followup
    if (x<end_study - eperiod) 0
    else if (x<=(end_study)) 1/(eperiod)
    else 0
  }
  CDF_x_followup<-function(x){ # CDF for x amount of followup
    if (x<end_study - eperiod) 0
    else if (x<=(end_study)) (x-(end_study - eperiod))/(eperiod)
    else 1
  }
  # control group
  toint<-function(u) sapply(u, function(t) (1-(1-CDF_x_followup(t))*survC$S(t))^(n/2-1)*P_x_followup(t)*survC$S(t))
  peval_c<-1-(n/2)*integrate(toint, lower = end_study - eperiod, upper = min(end_study, tau1))$value
  # trt group
  toint<-function(u) sapply(u, function(t) (1-(1-CDF_x_followup(t))*survT$S(t))^(n/2-1)*P_x_followup(t)*survT$S(t))
  peval_t<-1-(n/2)*integrate(toint, lower = end_study - eperiod, upper = min(end_study, tau1))$value
  return(peval_c*peval_t)
}
