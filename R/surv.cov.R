###############################
# surv.cov (one-arm) -- hidden
###############################

surv.cov <- function(df) {
  
  
  df$term <- ifelse((df$n.risk - df$n.event) == 0, 0, df$n.event / (df$n.risk * (df$n.risk - df$n.event)))
  
  # create empty covariance matrix of correct dimensions
  
  surv <- c(0,df$surv)
  term <- c(0, df$term)
  
  cov.mat = matrix(0, nrow = length(surv),length(surv))
  
  # populate control group covariance matrix
  
  for (i in 1:length(surv)) {
    for (j in 1:length(surv)) {
      t <- min(i, j)
      S.ti <- surv[i]
      S.tj <- surv[j]
      risk.sum <- sum(term[1:t])
      cov.mat[i,j] <- S.ti * S.tj * risk.sum
    }
  }
  
  Z = list()
  Z$cov <- cov.mat
  
  class(Z) = "surv.cov"
  return(Z)
  
  
}