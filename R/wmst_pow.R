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
#' @param n0 The total sample size in control group.
#' @param n1 The total sample size in treatment group.
#' @param power The desired power of the test.
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

wmst_pow <- function (survC, survT, censor = NA, eperiod, end_study, tau0, tau1, n0 = NA, n1 = NA, power = NA,
                      alpha = 0.05, two.sided = FALSE) {

  if (is.na(alpha))
    alpha <- ifelse(two.sided == T, 0.05, 0.025)
  if (!is.na(power) & (!is.na(n0) + !is.na(n1)) == 2)
    stop("One of n, power must be missing.")
  if (tau0 < 0)
    stop("Tau0 must be greater than 0.")
  if (tau1 <= tau0)
    stop("Tau1 must be greater than tau0.")
  if (is.na(n0) & is.na(n1)) {
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
    } else {
      find_n <- function(N) 1 - pnorm(qnorm(1 - alpha/2) - WMST_truediff/sqrt(evar(survT, eperiod, end_study, tau0, tau1, N) +
                                                             evar(survC, eperiod, end_study, tau0, tau1, N))) +
        pnorm(qnorm(alpha/2) - WMST_truediff/sqrt(evar(survT, eperiod, end_study, tau0, tau1, N) + evar(survC, eperiod, end_study, tau0, tau1, N))) - power
      if (find_n(10000) < 0)
        stop("Trial size would be more than 20,000 patients; please select different design parameters.")
    }
    N <- uniroot(find_n, lower = 1, upper = 10000)$root
    n0 <- ceiling(N)
    n1 <- ceiling(N)
  } else {
    NOTE <- paste("a total sample size of",n0 + n1)
    indicator <- "n0 + n1"

  }
  if (two.sided == F) {
    powerWMSTToverC <- powerWMST <- powfn(survC, survT,
                                          eperiod, end_study, tau0, tau1, n0, n1, alpha)
    powerWMSTCoverT <- NA
    powerLRToverC <- lrpow(survC, survT, eperiod,
                                 end_study, n0, n1, alpha, tau = end_study)
    powerLRCoverT <- NA
    powerLRtauToverC <- lrpow(survC, survT, eperiod,
                                    end_study, n0, n1, alpha, tau = tau1)
    powerLRtauCoverT <- NA

    lr.pow <- powerLRToverC
  } else {
    powerWMSTToverC <- powfn(survC, survT, eperiod, end_study,
                             tau0, tau1, n0, n1, alpha/2)
    powerWMSTCoverT <- powfn(survT, survC, eperiod, end_study,
                             tau0, tau1, n0, n1, alpha/2)
    powerWMST <- powerWMSTToverC + powerWMSTCoverT
    powerLRToverC <- lrpow(survC, survT, eperiod,
                                 end_study, n0, n1, alpha/2)
    powerLRCoverT <- lrpow(survT, survC, eperiod,
                                 end_study, n0 ,n1, alpha/2)
    powerLRtauToverC <- lrpow(survC, survT, eperiod,
                                    end_study, n0, n1, alpha/2, tau = tau1)
    powerLRtauCoverT <- lrpow(survT, survC, eperiod,
                                    end_study, n0, n1, alpha/2, tau = tau1)

    lr.pow <- powerLRToverC + powerLRCoverT
  }
  pKME <- wmst_eval(survC, survT, eperiod, end_study, tau0, tau1, n0, n1)
  to_ret <- list(NOTE = NOTE, indicator = indicator, tau0 = tau0, tau1 = tau1, two.sided = two.sided, n0 = n0, n1 = n1, powerWMST = powerWMST,
                 powerWMSTToverC = powerWMSTToverC,
                 powerWMSTCoverT = powerWMSTCoverT, powerLRToverC = powerLRToverC,
                 powerLRCoverT = powerLRCoverT, powerLRtauToverC = powerLRtauToverC,
                 powerLRtauCoverT = powerLRtauCoverT, pKME = pKME)

  
  Z = NULL
  Z$NOTE = NOTE
  Z$indicator = indicator
  Z$tau0 = tau0
  Z$tau1 = tau1
  Z$two.sided = two.sided
  Z$spec.pow = power
  Z$lr.pow = lr.pow
  Z$n0 = n0
  Z$n1 = n1
  Z$power = powerWMST
  Z$pkme = pKME
  class(Z) = "pow"
  return(Z)
}
