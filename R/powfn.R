################################
# powfn (two-arm) -- hidden
################################

powfn<-function(survC, survT, eperiod, end_study, tau0, tau1, n, alpha){
  WMST_truediff<-integrate(function(x) survT$S(x)-survC$S(x),
                           lower = tau0, upper = tau1)$value
  WMST_trueSE<-sqrt(evar(survT, eperiod, end_study, tau0, tau1, n/2)+
                      evar(survC, eperiod, end_study, tau0, tau1, n/2))
  return(1-pnorm(qnorm(1-alpha)-WMST_truediff/WMST_trueSE))
}