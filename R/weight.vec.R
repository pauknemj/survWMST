#################################
# weight.vec (one-arm) -- hidden
#################################


weight.vec <- function(df, tau0 = NULL, tau1 = NULL){
  
  #-- df -- dataframe of event and censoring times
  #-- tau0 -- lower truncation time
  #-- tau1 -- upper truncation time
  
  idx = which(df$time > tau0 & df$time < tau1)
  wk.time = sort(c(tau0, df$time[idx], tau1))
  time.diff <- diff(wk.time)
  if(tail(time.diff, n = 1) == 0) {
    time.diff <- head(time.diff,-1)
  }
  if(sum(idx) > 0){
    vid <- c(idx,tail(idx,1)+1)
  } else{vid = tail(which(df$time <= tau1),1) + 1}
  L <- rep(0, times = length(df$time) + 1)
  for(i in 1:length(vid)){
    L[vid[i]] <- time.diff[i]
  }
  
  Z = list()
  Z$weight <- as.matrix(L)
  
  class(Z) = "weight.vec"
  return(Z)
  
}