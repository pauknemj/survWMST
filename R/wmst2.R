#' @name wmst
#' @aliases wmst
#' @title Comparing window mean survival time
#' @description Performs two-sample comparisons using the window mean survival time (WMST) as a summary measure of the survival time distribution.
#' Three kinds of between-group contrast metrics (i.e., the difference in WMST, the ratio of WMST and the ratio of the window mean time lost (WMTL)) are computed.
#' The Greenwood plug-in estimator is used for the asymptotic variance. 
#' @usage  wmst(time, status, arm, tau0 = NULL, tau1 = NULL, alternative = "two.sided", alpha = 0.05)
#' @import survival
#' @param time The follow-up time for right censored data.
#' @param status The status indicator, 1=event, and 0=right censored.
#' @param arm The group indicator for comparison. The elements of this vector take either 1 or 0. Normally, 0=control group, 1=active treatment group.
#' @param tau0 A scaler value to specify the lower truncation time point of the window for the WMST calculation.
#' \code{tau0} needs to be at least 0 and smaller than tau1. When \code{tau0 = NULL}, the default value 0 is used (i.e., the lower bound for restricted mean survival analysis).
#' @param tau1 A scaler value to specify the upper truncation time point of the window for the WMST calculation.
#' \code{tau1} needs to be smaller than the minimum of the largest observed time in each of the two groups. When \code{tau1 = NULL}, the default value (i.e., the minimum of the largest observed time in each of the two groups) is used.
#' @param alternative Determines the sidedness of the alternative hypothesis. The default is "two.sided" which performs a two-tailed hypothesis test. If set to "greater" an upper-tailed test is performed. If set to "less" a lower-tailed test is performed.
#' @param alpha The default is 0.05. (1-\code{alpha}) confidence intervals are reported.
#' @details For more details, please see the package vignette: \code{browseVignettes(package = "survWMST")}
#' @return an object of class wmst.
#' @return \item{tau0}{the lower truncation time used in the analyses}
#' @return \item{tau1}{the upper truncation time used in the analyses}
#' @return \item{note}{a note regarding the truncation time}
#' @return \item{WMST.arm1}{WMST results in arm 1.}
#' @return \item{WMST.arm0}{WMST results in arm 0.}
#' @return \item{result}{Results of the WMST analyses.}
#' @references Note to self: figure out what to place here! (Definitely give Hajime Uno, Lu Tian, Miki Horiguchi, Angel Cronin, Chakib Battioui, James Bell credit)
#'
#' @author Mitchell Paukner



#'@export
#########################################
# wmst (2-arm) contrast (main function)
#########################################

wmst <- function (time, status, arm, tau0 = NULL, tau1 = NULL,
                   alternative = "two.sided", alpha = 0.05) {
  
  #==================================
  #  initial check
  #==================================
  
  #===== tau =====
  idx = arm == 0
  tt = time[idx]
  tau1.a.max = max(tt)
  idx = arm == 1
  tt = time[idx]
  tau1.b.max = max(tt)
  tau1_max = min(tau1.a.max, tau1.b.max)
  
  #--case 1: neither tau specified
  
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
  
  #--case 2: only tau0 unspecified
  
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
  
  #--case 3: only tau1 unspecified
  
  if (is.null(tau0) & !is.null(tau1)) {
    if (tau1 > tau1_max) {
      stop(paste("The truncation time, tau1, needs to be shorter than or equal to the minimum of the largest observed time of each of the two groups: ", 
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
  
  #--case 4: both tau specified
  
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
  
  #--ensuring alternative has been correctly specified
  
  if(alternative != "two.sided" & alternative != "less" & alternative != "greater"){
    stop(paste("The alternative must be : 'two.sided', 'less', or 'greater'."))    
  }
  
  Z = list()
  Z$tau0 = tau0
  Z$tau1 = tau1
  Z$note = NOTE
  if(tau0 != tau1){

      wk1 = wmst1(time[arm == 1], status[arm == 1], tau0, tau1, alpha)
      wk0 = wmst1(time[arm == 0], status[arm == 0], tau0, tau1, alpha)
      Z$WMST.arm1 = wk1
      Z$WMST.arm0 = wk0
      
      #--- contrast (WMST difference) ---
      
      wmst.diff.10 = wk1$wmst[1] - wk0$wmst[1]
      wmst.diff.10.se = sqrt(wk1$wmst.var + wk0$wmst.var)
      wmst.diff.10.low = wmst.diff.10 - qnorm(1 - alpha/2) * wmst.diff.10.se
      wmst.diff.10.upp = wmst.diff.10 + qnorm(1 - alpha/2) * wmst.diff.10.se
      if(alternative == "two.sided"){
        wmst.diff.pval = pnorm(-abs(wmst.diff.10)/wmst.diff.10.se) * 2
      } else if(alternative == "less"){
        wmst.diff.pval = pnorm(wmst.diff.10/wmst.diff.10.se)
      } else if(alternative == "greater"){wmst.diff.pval = pnorm(-wmst.diff.10/wmst.diff.10.se)}
      wmst.diff.result = c(wmst.diff.10, wmst.diff.10.low, 
                           wmst.diff.10.upp, wmst.diff.pval)
      
      #--- contrast (WMST ratio) ---
      
      wmst.log.ratio.10 = log(wk1$wmst[1]) - log(wk0$wmst[1])
      wmst.log.ratio.10.se = sqrt(wk1$wmst.var/wk1$wmst[1]/wk1$wmst[1] + 
                                    wk0$wmst.var/wk0$wmst[1]/wk0$wmst[1])
      wmst.log.ratio.10.low = wmst.log.ratio.10 - qnorm(1 - alpha/2) * wmst.log.ratio.10.se
      wmst.log.ratio.10.upp = wmst.log.ratio.10 + qnorm(1 - alpha/2) * wmst.log.ratio.10.se
      if(alternative == "two.sided"){
        wmst.log.ratio.pval = pnorm(-abs(wmst.log.ratio.10)/wmst.log.ratio.10.se) * 2
      } else if(alternative == "less"){
        wmst.log.ratio.pval = pnorm(wmst.log.ratio.10/wmst.log.ratio.10.se)
      } else if(alternative == "greater"){wmst.log.ratio.pval = pnorm(-wmst.log.ratio.10/wmst.log.ratio.10.se)}
      wmst.ratio.result = c(exp(wmst.log.ratio.10), exp(wmst.log.ratio.10.low), 
                            exp(wmst.log.ratio.10.upp), wmst.log.ratio.pval)
      
      #--- contrast (WMTL ratio  1/0) ---
      
      wmtl.log.ratio.10 = log(wk1$wmtl[1]) - log(wk0$wmtl[1])
      wmtl.log.ratio.10.se = sqrt(wk1$wmst.var/wk1$wmtl[1]/wk1$wmtl[1] + 
                                    wk0$wmst.var/wk0$wmtl[1]/wk0$wmtl[1])
      wmtl.log.ratio.10.low = wmtl.log.ratio.10 - qnorm(1 - alpha/2) * wmtl.log.ratio.10.se
      wmtl.log.ratio.10.upp = wmtl.log.ratio.10 + qnorm(1 - alpha/2) * wmtl.log.ratio.10.se
      if(alternative == "two.sided"){
        wmtl.log.ratio.pval = pnorm(-abs(wmtl.log.ratio.10)/wmtl.log.ratio.10.se) * 2
      } else if(alternative == "less"){
        wmtl.log.ratio.pval = pnorm(wmtl.log.ratio.10/wmtl.log.ratio.10.se)
      } else if(alternative == "greater"){wmtl.log.ratio.pval = pnorm(-wmtl.log.ratio.10/wmtl.log.ratio.10.se)}
      wmtl.ratio.result = c(exp(wmtl.log.ratio.10), exp(wmtl.log.ratio.10.low), 
                            exp(wmtl.log.ratio.10.upp), wmtl.log.ratio.pval)
      
      #--- results of WMST analysis ---
      
      out = rbind(wmst.diff.result, wmst.ratio.result, wmtl.ratio.result)
      if(tau0 != 0) {
        rownames(out) = c("WMST (arm=1)-(arm=0)",
                          "WMST (arm=1)/(arm=0)",
                          "WMTL (arm=1)/(arm=0)")
      } else if (tau0 == 0) {
        rownames(out) = c("RMST (arm=1)-(arm=0)",
                          "RMST (arm=1)/(arm=0)",
                          "RMTL (arm=1)/(arm=0)")
      }
      colnames(out) = c("Est.", paste("lower .", 
                                      round((1 - alpha) * 100, digits = 0), sep = ""), 
                        paste("upper .", round((1 - alpha) * 100, digits = 0), 
                              sep = ""), "p")
      
      #--- output of WMST analysis ---
      
      Z$result = out
      
      #--- WMST comparison by arm ---
      
      if (tau0 != 0){
        wmst = rbind(Z$WMST.arm1$wmst, Z$WMST.arm0$wmst)
        rownames(wmst) = c("WMST (arm=1)",
                           "WMST (arm=0)")
        Z$wmst = wmst
        
        #--- WMTL comparison by arm ---
        
        wmtl = rbind(Z$WMST.arm1$wmtl, Z$WMST.arm0$wmtl)
        rownames(wmtl) = c("WMTL (arm=1)",
                           "WMTL (arm=0)") 
        Z$wmtl = wmtl
      } else if (tau0 == 0){
        wmst = rbind(Z$WMST.arm1$wmst, Z$WMST.arm0$wmst)
        rownames(wmst) = c("RMST (arm=1)",
                           "RMST (arm=0)")
        Z$wmst = wmst
        
        #--- WMTL comparison by arm ---
        
        wmtl = rbind(Z$WMST.arm1$wmtl, Z$WMST.arm0$wmtl)
        rownames(wmtl) = c("RMTL (arm=1)",
                           "RMTL (arm=0)") 
        Z$wmtl = wmtl
      }
      

  } else if(tau0 == tau1){
    
    #--- Milestone Analysis ---
    
      wk1 = wmst1(time[arm == 1], status[arm == 1], tau0, tau1, alpha)
      wk0 = wmst1(time[arm == 0], status[arm == 0], tau0, tau1, alpha)
      Z$WMST.arm1 = wk1
      Z$WMST.arm0 = wk0
      
      #--- contrast (Milestone difference) ---
      
      wmst.diff.10 = wk1$wmst[1] - wk0$wmst[1]
      wmst.diff.10.se = sqrt(wk1$wmst.var + wk0$wmst.var)
      wmst.diff.10.low = wmst.diff.10 - qnorm(1 - alpha/2) * wmst.diff.10.se
      wmst.diff.10.upp = wmst.diff.10 + qnorm(1 - alpha/2) * wmst.diff.10.se
      if(alternative == "two.sided"){
        wmst.diff.pval = pnorm(-abs(wmst.diff.10)/wmst.diff.10.se) * 2
      } else if(alternative == "less"){
        wmst.diff.pval = pnorm(wmst.diff.10/wmst.diff.10.se)
      } else if(alternative == "greater"){wmst.diff.pval = pnorm(-wmst.diff.10/wmst.diff.10.se)}
      wmst.diff.result = c(wmst.diff.10, wmst.diff.10.low, 
                           wmst.diff.10.upp, wmst.diff.pval)
      
      #--- results of Milestone analysis ---
      
      out = rbind(wmst.diff.result)
      rownames(out) = c("Milestone (arm=1)-(arm=0)")
      colnames(out) = c("Est.", paste("lower .", 
                                      round((1 - alpha) * 100, digits = 0), sep = ""), 
                        paste("upper .", round((1 - alpha) * 100, digits = 0), 
                              sep = ""), "p")
      
      #--- output of Milestone analysis ---
      
      Z$result = out
      
      #--- WMST comparison by arm ---
      
      wmst = rbind(Z$WMST.arm1$wmst, Z$WMST.arm0$wmst)
      rownames(wmst) = c("Milestone (arm=1)",
                         "Milestone (arm=0)")
      colnames(wmst) = c("Est.", paste("lower .", 
                                       round((1 - alpha) * 100, digits = 0), sep = ""), 
                         paste("upper .", round((1 - alpha) * 100, digits = 0), 
                               sep = ""), "p")
      Z$wmst = wmst
      
      #--- WMTL comparison by arm ---
      
      wmtl = matrix(data=NA,nrow=2,ncol=4)
      rownames(wmtl) = c("WMTL (arm=1)",
                         "WMTL (arm=0)")
      colnames(wmtl) = c("Est.", paste("lower .", 
                                       round((1 - alpha) * 100, digits = 0), sep = ""), 
                         paste("upper .", round((1 - alpha) * 100, digits = 0), 
                               sep = ""), "p")
      Z$wmtl = wmtl
      
  }
  

  class(Z) = "wmst"
  Z
  
}

NULL