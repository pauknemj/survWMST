################################
# lrpow (two-arm) -- hidden
################################

lrpow<-function(survC, survT, eperiod = NA, end_study = NA, n0, n1, alpha, tau = NA){
  myk1<-round(eperiod, 3)
  myk2<-round(end_study - eperiod, 3)
  tot<-myk1+myk2
  
  L<-10000
  tis<-seq(from =0, to = tot, length.out = L)
  # Assess the hazard ratio and scaled hazard at the middle of each interval.
  
  xi<-survT$h(tis+tot/(2*(L-1)))/survC$h(tis+tot/(2*(L-1)))
  h_scl<-survC$h(tis+tot/(2*(L-1)))*tot/(L-1)
  
  nt<-nc<-rep(NA, L)
  nt[1]<-n0
  nc[1]<-n1

    for (i in 2:L){
      if (tis[i-1]<myk2){
        nc[i]<-nc[i-1]*(1-h_scl[i-1])
        nt[i]<-nt[i-1]*(1-h_scl[i-1]*xi[i-1])
      }
      else {
        nc[i]<-nc[i-1]*(1-h_scl[i-1]-(1/(L-(i-1))))
        nt[i]<-nt[i-1]*(1-h_scl[i-1]*xi[i-1]-(1/(L-(i-1))))
      }
    }

  nc[L]<-nt[L]<-0
  p<-nt/nc
  p[L]<-0
  D<-nc*h_scl + nt*h_scl*xi
  if (is.na(tau)) phi<-sum(D*(xi*p/(1+xi*p)-p/(1+p)))/sqrt(sum(D*p/(1+p)^2))
  else {
    if (tau<=0) stop('Tau must be greater than zero.')
    
    phi<-sum((D*(xi*p/(1+xi*p)-p/(1+p)))[tis<=tau])/sqrt(sum((D*p/(1+p)^2)[tis<=tau]))
  }
  
  return(pnorm(-qnorm(1-alpha)-phi))
}
