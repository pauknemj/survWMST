#################################
# window.cov (one-arm) -- hidden
#################################



window.cov <- function(time, status, tau = NULL){
  
  #-- time
  #-- status
  #-- tau -- dataframe of tau pairs for windows
  
  # Create data frames to import to surv.cov and weight.vec functions
  
  ft = survfit(Surv(time, status) ~ 1)
  df = with(ft, data.frame(time, n.risk, n.event, surv))
  if (sum(which(df$n.event == 0)) != 0) {
    df <- df[-which(df$n.event == 0), ]
  }
  
  
  
  # Get covariance matrix for each curve
  
  cov = surv.cov(df)$cov
  
  # create empty covariance matrix
  
  cov.mat = matrix(0, nrow = dim(tau)[1],dim(tau)[1])
  
  # fill covariance matrix
  
  
  for(i in 1:dim(cov.mat)[1]){
    for(j in 1:dim(cov.mat)[2]){
      L = weight.vec(df, tau0 = tau$tau0[i], tau1 = tau$tau1[i])$weight
      L.star = weight.vec(df, tau0 = tau$tau0[j], tau1 = tau$tau1[j])$weight
      cov.mat[i,j] = t(L) %*% cov %*% L.star
    }
  }  
  
  # correlation matrix
  
  cor.mat = cov2cor(cov.mat)
  
  # combine covariance matrices to get covariance matrix of difference
  
  
  Z = list()
  Z$cov = cov.mat
  Z$cor = cor.mat
  
  class(Z) = "window.cov"
  return(Z)
  
  
}