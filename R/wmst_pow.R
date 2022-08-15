#' @name wmst_pow
#' @aliases wmst_pow
#' @title Sample Size and Power for the Test of the Difference in Window Mean Survival Time
#' @description Determine the asymptotic power of the test of WMST under a given trial design, or calculate the samples size needed to achieve a desired power.
#' @usage  wmst_pow(survC, survT, censor = NA, eperiod, end_study, tau0, tau1, n = NA, power = NA, plot = FALSE, sim = F, M = 1000, method = "tau_star", alpha = 0.05, two.sided = FALSE)
#' @param survC the survival distribution of the control group, as a list in the form output by `wmst_surv`.
#' @param survT the survival distribution of the treatment group, as a list in the form output by `wmst_surv`.
#' @param censor The common CDF of the censoring distribution for treatment and control arms.
#' @param eperiod The duration of the enrollment period if subjects are assumed to accrue uniformly over the enrollment interval `(0,eperiod)`.
#' @param end_study The total duration of the study/follow-up period. 
#' @param tau0 A scaler vector to specify the lower truncation time points of the window for the WMST calculation.
#' \code{tau0} needs to be at least 0 and smaller than tau1. When \code{tau0 = NULL}, the default value 0 is used (i.e., the lower bound for restricted mean survival analysis).
#' @param tau1 A scaler value to specify the upper truncation time point of the window for the WMST calculation.
#' \code{tau1} needs to be smaller than the minimum of the largest observed time in each of the two groups. When \code{tau1 = NULL}, the default value (i.e., the minimum of the largest observed time in each of the two groups) is used.
#' @param n The total sample size for both groups. 1:1 randomization is assumed.
#' @param power The desired power of the test.
#' @param plot If TRUE, plots of the assumed survival distribution and power as a function of sample size are produced.
#' @param alpha The type I error level. Default is 0.05.
#' @param two.sided Determines whether the test is upper-tailed (`two.sided` = FALSE) or two-sided. The default is FALSE.
#' @details For more details, please see the package vignette: \code{browseVignettes(package = "survWMST")}
#' @return an object of class pow.
#' @return \item{n}{the user-specified n, or if n was left blank, the n needed to achieve the user-specified power.}
#' @return \item{powerWMST}{the user-specified power, or if power was left blank, the asymptotic power of the RMST test.
#' If `two.sided=FALSE`, `powerWMST` gives the power of rejecting in favor of the treatment arm.
#' If `two.sided=TRUE`, `powerWMST` gives the power of rejecting the null in favor of a two-sided alternative.}
#' @references Note to self: figure out what to place here! 
#'
#' @author Mitchell Paukner



#'@export
#########################################
# wmst_pow (2-arm) contrast (main function)
#########################################

wmst_pow <- function (survC, survT, censor = NA, eperiod, end_study, tau0, tau1, n = NA, power = NA, 
                     plot = FALSE, alpha = 0.05, two.sided = FALSE) {
  
  if (is.na(alpha)) 
    alpha <- ifelse(two.sided == T, 0.05, 0.025)
  if (is.na(power) + is.na(n) != 1) 
    stop("One of n, power must be missing.")
  if (tau0 < 0) 
    stop("Tau0 must be greater than 0.")
  if (tau1 <= tau0) 
    stop("Tau1 must be greater than tau0.")
  if (is.na(n)) {
    NOTE <- paste("a desired power of",power)
    indicator <- "power"
    WMST_truediff <- integrate(function(x) survT$S(x) - 
                                 survC$S(x), lower = tau0, upper = tau1)$value
    if (two.sided == F) {
      if (WMST_truediff < 0) 
        stop("True WMST in treatment arm is less than true WMST in control arm; cannot design a trial to show control is superior.")
      WMST_trueSE <- WMST_truediff/(qnorm(1 - alpha) - 
                                      qnorm(1 - power))
      find_n <- function(N) sqrt(evar(survT, eperiod, end_study, 
                                      tau0, tau1, N) + evar(survC, eperiod, end_study, tau0, tau1, N)) - WMST_trueSE
      if (find_n(10000) > 0) 
        stop("Trial size would be more than 20,000 patients; please select different design parameters.")
    }
    else {
      find_n <- function(N) 1 - pnorm(qnorm(1 - alpha/2) - WMST_truediff/sqrt(evar(survT, eperiod, end_study, tau0, tau1, N) + 
                                                             evar(survC, eperiod, end_study, tau0, tau1, N))) + 
        pnorm(qnorm(alpha/2) - WMST_truediff/sqrt(evar(survT, eperiod, end_study, tau0, tau1, N) + evar(survC, eperiod, end_study, tau0, tau1, N))) - power
      if (find_n(10000) < 0) 
        stop("Trial size would be more than 20,000 patients; please select different design parameters.")
    }
    N <- uniroot(find_n, lower = 1, upper = 10000)$root
    n <- 2 * ceiling(N)
  }
  else {
    NOTE <- paste("a total sample size of",n)
    indicator <- "n"
    if (n%%2 != 0) 
      n <- ifelse(floor(n)%%2 == 0, floor(n), floor(n) - 
                    1)
  }
  if (two.sided == F) {
    powerWMSTToverC <- powerWMST <- powfn(survC, survT, 
                                          eperiod, end_study, tau0, tau1, n, alpha)
    powerWMSTCoverT <- NA
    powerLRToverC <- lrpow(survC, survT, eperiod, 
                                 end_study, n, alpha, tau = end_study)
    powerLRCoverT <- NA
    powerLRtauToverC <- lrpow(survC, survT, eperiod, 
                                    end_study, n, alpha, tau = tau1)
    powerLRtauCoverT <- NA
    
    lr.pow <- powerLRToverC
  }
  else {
    powerWMSTToverC <- powfn(survC, survT, eperiod, end_study, 
                             tau0, tau1, n, alpha/2)
    powerWMSTCoverT <- powfn(survT, survC, eperiod, end_study, 
                             tau0, tau1, n, alpha/2)
    powerWMST <- powerWMSTToverC + powerWMSTCoverT
    powerLRToverC <- lrpow(survC, survT, eperiod, 
                                 end_study, n, alpha/2)
    powerLRCoverT <- lrpow(survT, survC, eperiod, 
                                 end_study, n, alpha/2)
    powerLRtauToverC <- lrpow(survC, survT, eperiod, 
                                    end_study, n, alpha/2, tau = tau1)
    powerLRtauCoverT <- lrpow(survT, survC, eperiod, 
                                    end_study, n, alpha/2, tau = tau1)
    
    lr.pow <- powerLRToverC + powerLRCoverT
  }
  pKME <- wmst_eval(survC, survT, eperiod, end_study, tau0, tau1, n)
  to_ret <- list(NOTE = NOTE, indicator = indicator, tau0 = tau0, tau1 = tau1, two.sided = two.sided, n = n, powerWMST = powerWMST, 
                 powerWMSTToverC = powerWMSTToverC, 
                 powerWMSTCoverT = powerWMSTCoverT, powerLRToverC = powerLRToverC, 
                 powerLRCoverT = powerLRCoverT, powerLRtauToverC = powerLRtauToverC, 
                 powerLRtauCoverT = powerLRtauCoverT, pKME = pKME)
  
  #plot 
  
  if (plot == T) {
    
    # par(mfrow = c(2, 2), mar = c(5, 4, 1, 2))
    
    par(mfrow = c(1, 2))
    
    n_vec <- seq(from = n * 0.5, to = n * 1.5, length.out = 20)
    
    if (eperiod == 0){
      eperiod_vec <- seq(from = 0, to = end_study - eperiod, length.out = 20) 
    } else{eperiod_vec <- seq(from = max(0, tau1- (end_study - eperiod)), to = eperiod * 1.5, length.out = 20)}
    
    if (end_study - eperiod == 0){
      k2_vec <- seq(from = 0, to = eperiod, length.out = 20)
    } else{k2_vec <- seq(from = max(0, tau1- eperiod), to = (end_study - eperiod) * 1.5, length.out = 20)}
    
    tau_vec <- seq(from = tau1* 0.5, to = min(end_study, tau1*1.5), length.out = 20)
    
    plotsurvdef(survC, survT, xupper = end_study)
    
    legend((end_study) * 0.5, 0.25, c("Control", "Treatment"), 
           cex = 0.5, lty = c(1, 2), bty = "n")
    
    if (two.sided == F) {
      plot(n_vec, sapply(n_vec, function(x) powfn(survC, 
                                                  survT, eperiod, end_study, tau0, tau1, x, alpha)), xlab = "Sample size", 
           ylab = "Power", ylim = c(0, 1), type = "l")
      
      lines(n_vec, sapply(n_vec, function(x) powfn(survC, 
                                                   survT, eperiod, end_study, tau0 = 0, tau1, x, alpha)), col = "red")
      
      lines(n_vec, sapply(n_vec, function(x) lrpow(survC, 
                                                         survT, eperiod, end_study, x, alpha)), col = "blue")
      
      
      legend(n_vec[10], 0.25, c("WMST", "RMST", 
                                "Log-rank"), col = c(1, 2, 4), cex = 0.5, 
             lty = 1, bty = "n")

    }
    else {
      plot(n_vec, sapply(n_vec, function(x) powfn(survC, 
                                                  survT, eperiod, end_study, tau0, tau1, x, alpha/2)), xlab = "Sample size", 
           ylab = "Power", ylim = c(0, 1), type = "l")
      
      lines(n_vec, sapply(n_vec, function(x) powfn(survT, 
                                                   survC, eperiod, end_study, tau0, tau1, x, alpha/2)), lty = 2)
      
      lines(n_vec, sapply(n_vec, function(x) powfn(survC, 
                                                   survT, eperiod, end_study, tau0 = 0, tau1, x, alpha/2)), col = "red")
      
      lines(n_vec, sapply(n_vec, function(x) powfn(survT, 
                                                   survC, eperiod, end_study, tau0 = 0, tau1, x, alpha/2)), col = "red", 
            lty = 2)
      
      lines(n_vec, sapply(n_vec, function(x) lrpow(survC, 
                                                         survT, eperiod, end_study, x, alpha/2)), col = "blue")
      
      lines(n_vec, sapply(n_vec, function(x) lrpow(survT, 
                                                         survC, eperiod, end_study, x, alpha/2)), col = "blue", 
            lty = 2)
  
      
      legend(n_vec[10], 0.25, c("WMST", "RMST", 
                                "Log-rank"), col = c(1, 2, 4), cex = 0.5, 
             lty = 1, bty = "n")
      
    }
  }
  
  Z$NOTE = NOTE
  Z$indicator = indicator
  Z$tau0 = tau0
  Z$tau1 = tau1
  Z$two.sided = two.sided
  Z$spec.pow = power
  Z$lr.pow = lr.pow
  Z$n = n
  Z$power = powerWMST
  Z$pkme = pKME
  class(Z) = "pow"
  return(Z)
}
