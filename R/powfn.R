################################
# powfn (two-arm) -- hidden
################################

powfn<-function(survC, survT, eperiod, end_study, tau0, tau1, n0, n1, alpha){
  WMST_truediff<-integrate(function(x) survT$S(x)-survC$S(x),
                           lower = tau0, upper = tau1)$value
  WMST_trueSE<-sqrt(evar(survT, eperiod, end_study, tau0, tau1, n0)+
                      evar(survC, eperiod, end_study, tau0, tau1, n1))
  return(1-pnorm(qnorm(1-alpha)-WMST_truediff/WMST_trueSE))
}
