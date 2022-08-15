################################
# evar (one-arm) -- hidden
################################

# evar<-function(survo, censor = NA,  eperiod, end_study, tau0, tau1, n){
#
#     P_censor <- function(t,x) {
#       if (t<=x) 1
#       else if (x<t & t<=(eperiod+x))  1 - (t-x)/eperiod
#       else if (eperiod+x <t & t<= end_study) 0
#       else if (eperiod+x <t & all.equal(t, end_study)) 0}
#
#   toint<-function(u) sapply(u, function(x) integrate(survo$S, lower = max(tau0,x), upper = tau1)$value^2*survo$h(x)/(n*survo$S(x)*(1 - P_censor(end_study, x))))
#   integrate(toint, lower = 0, upper = tau1)$value
# }

evar<-function(survdefo,  eperiod, end_study, tau0, tau1, n){
  P_censor<-function(t,x) {
    if (t<=x) 1
    else if (x<t & t<=(eperiod+x))  1 - (t-x)/eperiod
    else if (eperiod+x <t & t<= end_study) 0
    else if (eperiod+x <t & all.equal(t, end_study)) 0}
  toint<-function(u) sapply(u, function(x) integrate(survdefo$S, lower = max(tau0,x), upper = tau1)$value^2*survdefo$h(x)/(n*survdefo$S(x)*(1 - P_censor(end_study, x))))
  integrate(toint, lower = 0, upper = tau1)$value
}

# evar<-function(survdefo,  k1, k2, tau0, tau1, n){
#   P_censor<-function(t,x) {
#     if (t<=x) 1
#     else if (x<t & t<=(k1+x))  1 - (t-x)/k1
#     else if (k1+x <t & t<= k1+k2) 0
#     else if (k1+x <t & all.equal(t, k1+k2)) 0}
#   toint<-function(u) sapply(u, function(x) integrate(survdefo$S, lower = max(tau0,x), upper = tau1)$value^2*survdefo$h(x)/(n*survdefo$S(x)*(1 - P_censor(k1 + k2, x))))
#   integrate(toint, lower = 0, upper = tau1)$value
# }
