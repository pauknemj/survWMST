#' @name combo
#' @aliases wmst2
#' @import mvtnorm
#' @title Comparing window mean survival time
#' @description Performs two-sample comparisons using the window mean survival time (WMST) as a summary measure of the survival time distribution.
#' Three kinds of between-group contrast metrics (i.e., the difference in WMST, the ratio of WMST and the ratio of the window mean time lost (WMTL)) are computed.
#' The Greenwood plug-in estimator is used for the asymptotic variance.
#' @usage  combo(time, status, arm, type = NULL, combination = FALSE, round = FALSE, method = c("maximum", "average", "chi.sq")[1], tau0 = NULL, tau1 = NULL, alternative = c("two.sided", "greater", "less")[1], alpha = 0.05)
#' @param time The follow-up time for right censored data.
#' @param status The status indicator, 1=event, and 0=right censored.
#' @param arm The group indicator for comparison. The elements of this vector take either 1 or 0. Normally, 0=control group, 1=active treatment group.
#' @param type The type of window selection method to be applied.
#' @param combination Determines whether the windows should be selected as every combination (TRUE) of the choices of tau0 and tau1 or not (FALSE).
#' @param round Determines whether the choices for tau0 and tau1 will be rounded to the nearest whole number (TRUE) or will remain as chosen (FALSE).
#' @param method Selects how the test statistics from each window will be combined.
#' @param tau0 A scaler vector to specify the lower truncation time points of the window for the WMST calculation.
#' \code{tau0} needs to be at least 0 and smaller than tau1. When \code{tau0 = NULL}, the default value 0 is used (i.e., the lower bound for restricted mean survival analysis).
#' @param tau1 A scaler vector to specify the upper truncation time points of the window for the WMST calculation.
#' \code{tau1} needs to be smaller than the minimum of the largest observed time in each of the two groups. When \code{tau1 = NULL}, the default value (i.e., the minimum of the largest observed time in each of the two groups) is used.
#' @param alternative Determines the sidedness of the alternative hypothesis. The default is "two.sided" which performs a two-tailed hypothesis test. If set to "greater" an upper-tailed test is performed. If set to "less" a lower-tailed test is performed.
#' @param alpha The default is 0.05. (1-\code{alpha}) confidence intervals are reported.
#' @details For more details, please see the package vignette: \code{browseVignettes(package = "survWMST")}
#' @return an object of class combo.
#' @return \item{pvalue}{the pvalue produced from the versatile test of multiple windows.}
#' @return \item{tau}{the combination of various choices of tau0 and tau1 for windows.}
#' @references Note to self: figure out what to place here!
#'
#' @author Mitchell Paukner



#'@export
#########################################
# combo (2-arm) contrast (main function)
#########################################

combo <- function(time, status, arm, type = NULL, combination = FALSE, round = FALSE, method = c("maximum", "average", "chi.sq")[1],
                  tau0 = NULL, tau1 = NULL, alternative = c("two.sided", "greater", "less")[1], alpha = 0.05){


  idx = arm == 0
  tt = time[idx]
  tau1.a.max = max(tt)
  idx = arm == 1
  tt = time[idx]
  tau1.b.max = max(tt)
  tau1_max = min(tau1.a.max, tau1.b.max)


  # Get choices of tau0 and tau1

  if (is.null(tau0) & is.null(tau1)){
    tau = tau.pairs(time, status, type, round, combination, tau1_max = tau1_max)

    if(type == "quarter"){
      if(combination){
        NOTE = paste("The time quarters method was used for window selection \n with", tau1_max, "as an upper bound (with combination).")
      } else{
        NOTE = paste("The time quarters method was used for window selection \n with", tau1_max, "as an upper bound (without combination).")
      }
    } else if(type == "quantile"){
      if(combination){
        NOTE = paste("The time quantiles method was used for window selection \n with", tau1_max, "as an upper bound (with combination).")
      } else{
        NOTE = paste("The time quantiles method was used for window selection \n with", tau1_max, "as an upper bound (without combination).")
      }
    } else if(type == "quartile"){
      if(combination){
        NOTE = paste("The time event time quartiles method was used for window selection \n with", tau1_max, "as an upper bound (with combination).")
      } else{
        NOTE = paste("The event time quartiles method was used for window selection \n with", tau1_max, "as an upper bound (without combination).")
      }
    } else if(type == "survquart"){
      if(combination){
        NOTE = paste("The pseudo survival quartiles method was used for window selection \n with", tau1_max, "as an upper bound (with combination).")
      } else{
        NOTE = paste("The pseudo survival quartiles method was used for window selection \n with", tau1_max, "as an upper bound (without combination).")
      }
    }

  } else{
    tau = data.frame("tau0" = tau0, "tau1" = tau1)
    if(length(tau0) == length(tau1)){
      NOTE = paste("The windows were specified by the user.")
    } else{
      stop(paste("The vectors for tau0 and tau1 must be of equal length."))
    }
  }




  # Covariance matrix for both arms

  cov0 = window.cov(time[arm == 0], status[arm == 0], tau)$cov
  cov1 = window.cov(time[arm == 1], status[arm == 1], tau)$cov


  # Covariance of difference

  cov.diff = cov0 + cov1
  cor.diff = cov2cor(cov.diff)

  # get difference in window means
  Z = list()
  diff = NULL
  se = NULL
  test.stat = NULL

  for(i in 1:dim(tau)[1]){

    diff[i] = window.mean(time[arm == 1], status[arm == 1], tau0 = tau$tau0[i], tau1 = tau$tau1[i])$wmst -
      window.mean(time[arm == 0], status[arm == 0], tau0 = tau$tau0[i], tau1 = tau$tau1[i])$wmst
    se[i] = sqrt(cov.diff[i,i])
    test.stat[i] = diff[i]/sqrt(cov.diff[i,i])


  }



  if(method == "maximum"){
    if(alternative == "two.sided"){
      max.stat = max(abs(test.stat))
      upper = 1 - pmvnorm(lower = rep(-Inf, dim(cor.diff)[1]), upper = rep(max.stat,dim(cor.diff)[1]), corr = cor.diff)[1]
      lower = pmvnorm(lower = rep(-Inf, dim(cor.diff)[1]), upper = rep(max.stat,dim(cor.diff)[1]), corr = cor.diff)[1]
      pvalue = 2*min(upper,lower)
      max.window.tau0 = tau$tau0[which(abs(test.stat) == max.stat)]
      max.window.tau1 = tau$tau1[which(abs(test.stat) == max.stat)]
    } else if(alternative == "less"){
      max.stat = min(test.stat)
      pvalue = pmvnorm(lower = rep(-Inf, dim(cor.diff)[1]), upper = rep(max.stat,dim(cor.diff)[1]), corr = cor.diff)[1]
      max.window.tau0 = tau$tau0[which(test.stat == max.stat)]
      max.window.tau1 = tau$tau1[which(test.stat == max.stat)]
    } else if(alternative == "greater"){
      max.stat = max(test.stat)
      pvalue = 1 - pmvnorm(lower = rep(-Inf, dim(cor.diff)[1]), upper = rep(max.stat,dim(cor.diff)[1]), corr = cor.diff)[1]
      max.window.tau0 = tau$tau0[which(test.stat == max.stat)]
      max.window.tau1 = tau$tau1[which(test.stat == max.stat)]
    }
  } else if(method == "chi.sq"){
    chi.sq = t(test.stat) %*% solve(cor.diff) %*% test.stat
    pvalue = pchisq(chi.sq, df = dim(cor.diff)[1], lower.tail = FALSE)
  } else if(method == "average"){
    avg.stat = mean(test.stat)
    avg.var = sum(cor.diff)/dim(cor.diff)[1]^2
    z.score = avg.stat/sqrt(avg.var)
    if(alternative == "two.sided"){
      pvalue = pnorm(-abs(z.score))*2
    } else if(alternative == "less"){
      pvalue = pnorm(z.score)
    } else if(alternative == "greater"){
      pvalue = pnorm(-z.score)
    }
  }

  out = cbind(tau$tau0, tau$tau1, diff, se, test.stat)

  colnames(out) = c("Tau0", "Tau1", "Est.", "se", "Z-score")


  Z$result = out
  Z$alternative = alternative
  Z$method = method
  Z$NOTE = NOTE
  Z$p = pvalue
  Z$tau = tau
  Z$diff = diff
  class(Z) = "combo"
  return(Z)


}
