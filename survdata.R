#' @name survdata
#' @aliases survdata
#' @title Simulate Exponential Survival Data
#' @description Creates a data frame of survival data based on exponential distributions and Bernoulli censoring for a 2 arm trial. 
#' @usage  rmst2(n, lambda1, lambda2, lambda3, censor = 0, piecewise = FALSE, diverge_time = NULL, end_study = NULL)
#' @param n The sample size in each trial arm.
#' @param lambda1 The rate of the exponential distribution in the first arm.
#' @param lambda2 The rate of the exponential distribution in the second arm if piecewise = FALSE, or the rate of the exponential distribution of the second arm before diverge_time if piecewise = TRUE.
#' @param lambda3 the rate of the exponential distribution of the second arm after diverge_time if piecewise = TRUE.
#' @param censor The probability a sample time will be deemed a censored observation.
#' \code{censor} needs to be between 0 and 1.
#' @param piecewise Determines if the second trial arm follows a piecewise exponential distribution.
#' @param diverge_time Determines the point in time in which the piecewise exponential curve changes rates from lambda 2 to lambda 3.
#' @param end_study Creates a point at which all the simulated times are censored after.
#' @return dataframe with oberved times, censoring indicators (0 for censoring, 1 for events), and trial arm.
#' @references Note to self: figure out what to place here! (Definitely give Hajime Uno, Lu Tian, Miki Horiguchi, Angel Cronin, Chakib Battioui, James Bell credit)
#'
#' @author Mitchell Paukner



#'@export
#########################################
# survdata (function)
#########################################


survdata <- function(n, lambda1, lambda2, lambda3, censor = 0, piecewise = FALSE,
                          diverge_time = NULL, end_study = NULL){
  
  if(piecewise){
    curve1 <- rexp(n, rate = lambda1)
    curve2 <- rpwexp(n, rate = c(lambda2,lambda3), time = c(0,diverge_time))
  } else {
    curve1 <- rexp(n, rate = lambda1)
    curve2 <- rexp(n, rate = lambda2)
  }
  arm <- c(rep(c("0","1"), each = n))
  time <- c(curve1, curve2)
  
  if(censor < 1){
    event <- NULL
    
    for(j in 1:(2*n)){
      if(!is.null(end_study)){
        if(time[j] <= end_study){event[j] = rbinom(1, 1, 1 - censor)}
        else{event[j] = 0; time[j] = end_study}
      } else{event[j] = 1}
    }
  } else{
    event <- c(rep(1, times = 2 * n))
  }
  
  return(data.frame("time" = time, "event" = event, "arm" = arm))
  
}
