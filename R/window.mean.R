##################################
# window.mean (one-arm) -- hidden
##################################

window.mean <- function(time, status, tau0 = NULL, tau1 = NULL){
  
  #-- time
  #-- status
  #-- tau0 -- lower truncation time
  #-- tau1 -- upper truncation time  
  
  tau1_max <- max(time)
  if (!is.null(tau0) & !is.null(tau1)) {
    if (tau0 < 0) {
      stop(paste("The truncation time, tau0, needs to be larger than or equal to 0."))
    }     
    if (tau1 < tau0) {
      stop(paste("The truncation time, tau1, needs to be larger than or equal to tau0: ", 
                 round(tau0, digits = 3)))
    }  
    if (tau1 > tau1_max) {
      stop(paste("The truncation time, tau1, needs to be shorter than or equal to the minimum of the largest observed time on each of the two groups: ", 
                 round(tau1_max, digits = 3)))
    }
    if (tau1 >= tau0 & tau1 <= tau1_max) {
      NOTE = paste("The truncation times: tau0 =", 
                   tau0, " and tau1 =", tau1, " were specified.")
    }    
  }
  if (!is.null(tau0) & is.null(tau1)) {
    if (tau0 < 0) {
      stop(paste("The truncation time, tau0, needs to be larger than or equal to 0."))
    }    
    tau1 = tau1_max
    if (tau1 < tau0) {
      stop(paste("The truncation time, tau1, needs to be larger than or equal to tau0: ", 
                 round(tau0, digits = 3)))
    }    
    if (tau1 >= tau0) {
      NOTE = (paste("The truncation time, tau1, was not specified. Thus, the default tau1 (the minimum of the largest observed time on each of the two groups)", 
                    round(tau_max, digits = 3), " is used, while tau0 was specified as: tau0 =", round(tau0, digits = 3)))
    }
  }
  if (is.null(tau0) & !is.null(tau1)) {
    if (tau1 > tau1_max) {
      stop(paste("The truncation time, tau1, needs to be shorter than or equal to the minimum of the largest observed time on each of the two groups: ", 
                 round(tau1_max, digits = 3)))
    }    
    tau0 = 0
    if (tau1 < tau0) {
      stop(paste("The truncation time, tau1, needs to be larger than or equal to tau0: ", 
                 round(tau0, digits = 3)))
    }    
    if (tau1 >= tau0) {
      NOTE = (paste("The truncation time, tau0, was not specified. Thus, the default tau0 of 0 is used, while tau1 was specified as: tau1 =", round(tau1, digits = 3)))
    }
  }  
  if (is.null(tau0) & is.null(tau1)) {
    tau0 = 0
    tau1 <- tau1_max
    if (tau1 < tau0) {
      stop(paste("The truncation time, tau1, needs to be larger than or equal to tau0: ", 
                 round(tau0, digits = 3)))
    }    
    if (tau1 >= tau0) {
      NOTE = (paste("The truncation times, tau0 and tau1, were not specified. Thus, the default tau0 = 0 and the default tau1 (the minimum of the largest observed time on each of the two groups)", round(tau1_max, digits = 3)," are used."))
    }
  } 
  
  
  # compute window mean
  
  ft = survfit(Surv(time, status) ~ 1)
  df <- with(ft, data.frame(time, n.risk, n.event, surv))
  if (sum(df$n.event) != length(df$n.event)) {
    if(sum(which(df$n.event == 0)) != 0){
      df <- df[-which(df$n.event == 0), ]
    }
  }
  if(tau0 != tau1){
    idx = which(df$time >= tau0 & df$time <= tau1)
  } else {idx = tail(which(df$time <= tau1),1)}
  wk.time = sort(c(tau0, df$time[idx], tau1))
  wk.surv = df$surv[idx]
  wk.n.risk = df$n.risk[idx]
  wk.n.event = df$n.event[idx]
  time.diff <- diff(wk.time)
  if(min(df$time) < tau0){
    if(sum(idx) == 0){
      areas <- time.diff * df$surv[which.max(df$time[which(df$time <= tau0)])]
    } else{areas <- time.diff * c(df$surv[min(idx) - 1], wk.surv)}
  } else{areas <- time.diff * c(1, wk.surv)}
  if(tau0 < tau1){
    wmst = sum(areas)
  } else if(tau0 == tau1){
    wmst = wk.surv
  }
  
  Z = list()
  Z$note = NOTE
  Z$wmst = wmst
  
  class(Z) = "window.mean"
  return(Z)
  
}

