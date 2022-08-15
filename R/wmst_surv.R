#' Create a Survdef Object for a Piecewise Exponential Distribution
#'
#' Create a new object which stores user-specified survival distribution
#' information in the format needed for the main function, `RMSTpow`. `survdef`
#' is used when the user wishes to specify a piecewise expoenential survival
#' distribution. Either the hazard on fixed intervals or survival probabilities
#' at fixed times can be specified.
#'
#' @param haz a vector of hazards of length<=10. If a single hazard is
#' specified, the survival distribution is exponential with the specified
#' hazard. If `haz` has length>2, the survival distribution has constant hazard
#' equal to `haz` over the intervals `[0, t_1), [t_1, t_2), \dots, [t_m, Inf)`
#' where `t_i` are the entries of `times` and `times` has length `m`. One of
#' `haz`, `surv` must be specified.
#' @param surv a vector of survival probabilities of length<=10 corresponding to
#' `times`. If `surv` is specified, the survival distribution has constant hazard
#' over the intervals `[0, t_1), [t_1, t_2), \dots, [t_m, Inf)`
#' where `t_i` are the entries of `times` and `times` has length `m+1`. The hazards
#' are calculated so that the curve passes through each entry in `surv` at the
#' corresponding time from `times`. One of `haz`, `surv` must be specified.
#' @param times a vector of the same length as surv (if surv is specified) or
#' one element shorter than haz (if haz is specified). No times term is required
#'  if a single hazard is specified in haz.
#'
#' @return a list with components:
#' \item{S}{a vectorized function that takes time as input and returns the survival probability at that time}
#' \item{h}{a vectorized function that takes time as input and returns the hazard at that time}
#' @export

################################
# wmst_surv (one-arm)
################################

wmst_surv <- function (haz = NA, surv = NA, times = NA) {

  if (is.na(haz[1]) + is.na(surv[1]) != 1)
    stop("You must specify either haz or surv.")
  nti <- length(times)
  if (!is.na(surv[1])) {
    if (length(surv) != nti | is.na(times[1]))
      stop("Length of times and surv must match.")
    if (nti == 1)
      h <- function(x) sapply(x, function(y) log(1/surv)/times)
    else {
      if (nti > 10)
        stop("You can only specify up to 10 survival probabilities.")
      times10 <- c(times[-nti], rep(Inf, 10 - nti))
      myh <- log(c(1, surv[-nti])/c(surv))/diff(c(0, times))
      h <- function(x) {
        sapply(x, function(t) {
          if (0 <= t & t < times10[1])
            myh[1]
          else if (times10[1] <= t & t < times10[2])
            myh[2]
          else if (times10[2] <= t & t < times10[3])
            myh[3]
          else if (times10[3] <= t & t < times10[4])
            myh[4]
          else if (times10[4] <= t & t < times10[5])
            myh[5]
          else if (times10[5] <= t & t < times10[6])
            myh[6]
          else if (times10[6] <= t & t < times10[7])
            myh[7]
          else if (times10[7] <= t & t < times10[8])
            myh[8]
          else if (times10[8] <= t & t < times10[9])
            myh[9]
          else if (times10[9] <= t)
            myh[10]
        })
      }
    }
  }
  else if (!is.na(haz[1])) {
    if (is.na(times[1]) & length(haz) != 1)
      stop("Without times specified, the only valid input is a single hazard.")
    if (length(haz) == 1)
      h <- function(x) sapply(x, function(y) haz)
    else {
      if (length(haz) > 10)
        stop("You can only specify up to 10 hazards.")
      if (length(haz) != nti + 1)
        stop("Length of haz must be one longer than times.")
      times10 <- c(times, rep(Inf, 9 - nti))
      h <- function(x) {
        sapply(x, function(t) {
          if (0 <= t & t < times10[1])
            haz[1]
          else if (times10[1] <= t & t < times10[2])
            haz[2]
          else if (times10[2] <= t & t < times10[3])
            haz[3]
          else if (times10[3] <= t & t < times10[4])
            haz[4]
          else if (times10[4] <= t & t < times10[5])
            haz[5]
          else if (times10[5] <= t & t < times10[6])
            haz[6]
          else if (times10[6] <= t & t < times10[7])
            haz[7]
          else if (times10[7] <= t & t < times10[8])
            haz[8]
          else if (times10[8] <= t & t < times10[9])
            haz[9]
          else if (times10[9] <= t)
            haz[10]
        })
      }
    }
  }
  if (is.na(times[1]))
    S <- function(t) sapply(t, function(z) exp(-haz * z))
  else S <- function(t) sapply(t, function(z) {
    width <- diff(c(0, times[times < z], z))
    H <- sum(width * h(c(0, times[times < z])))
    return(exp(-H))
  })
  return(list(S = S, h = h))
}
