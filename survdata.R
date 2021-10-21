#' @name survdata
#' @aliases survdata
#' @title Simulate Exponential Survival Data
#' @description Creates a data frame of survival data based on exponential distributions and Bernoulli censoring for a 2 arm trial. 
#' @usage  survdata <- function(n, lambda1, lambda2, censor = 0, piecewise = FALSE, diverge_time = NULL, erate = NULL, eperiod = NULL, end_study = NULL)
#' @param n The sample size in each trial arm.
#' @param lambda1 The rate of the exponential distribution in the first arm.
#' @param lambda2 The rate of the exponential distribution in the second arm if piecewise = FALSE, or the rates of the piecewise exponential distribution of the second arm before diverge_time if piecewise = TRUE.
#' @param censor The probability a sample time will be deemed a censored observation.
#' \code{censor} needs to be between 0 and 1.
#' @param piecewise Determines if the second trial arm follows a piecewise exponential distribution.
#' @param diverge_time Determines the point in time in which the piecewise exponential curve changes rates.
#' @param erate Vector of enrollment rates following piecewise exponential distribution.
#' @param eperiod Vector of length +1 to erate which determines the times at which the enrollment period begins (0), when the erates change, and when the enrollment period ends.
#' @param end_study Creates a point at which all the simulated times are censored after.
#' @return dataframe with observed times, calendar times, enrollment times, censoring indicators (0 for censoring, 1 for events), and trial arm.
#'
#' @author Mitchell Paukner



#'@export
#########################################
# survdata (function)
#########################################


survdata <- function(n, lambda1, lambda2, censor = 0, piecewise = FALSE, diverge_time = NULL, 
                     erate = NULL, eperiod = NULL, end_study = NULL){
  
  if(!is.null(erate) & !is.null(eperiod)){
    etime <- rpwexp(2*n, rate = c(erate), time = head(eperiod, -1)) # enrollment time 
    i <- 1
    while(i <= 2*n){
      if(etime[i] <= tail(eperiod,1)){
        i <- i + 1
      } else{
        etime[i] <- rpwexp(1, rate = c(erate), time = head(eperiod, -1))
      }
    }
  } else{
    etime <- c(rep(0, times = 2*n))
  }
  
  if(piecewise){
    survt1 <- rexp(n, rate = lambda1) # survival time control
    survt2 <- rpwexp(n, rate = c(lambda2), time = c(0,diverge_time)) # survival time trt
    calt1 <- survt1 + etime[1:n] # calendar time control
    calt2 <- survt2 + etime[(n+1):(2*n)] # calendar time trt
  } else {
    survt1 <- rexp(n, rate = lambda1)
    survt2 <- rexp(n, rate = lambda2[1])
    calt1 <- survt1 + etime[1:n]
    calt2 <- survt2 + etime[(n+1):(2*n)]
  }
  arm <- c(rep(c("0","1"), each = n))
  time <- c(survt1, survt2)
  calt <- c(calt1, calt2)
  
  if(censor < 1){
    event <- NULL
    
    for(j in 1:(2*n)){
      if(!is.null(end_study)){
        if(calt[j] <= end_study){event[j] = rbinom(1, 1, 1 - censor)}
        else{event[j] = 0; time[j] <- ifelse(time[j] <= end_study, time[j], end_study)}
      } else{event[j] = rbinom(1, 1, 1 - censor)}
    }
  } else{
    event <- c(rep(1, times = 2 * n))
  }
  
  return(data.frame("time" = time, "calendar" = calt, "enroll" = etime, "event" = event, "arm" = arm))
  
}

