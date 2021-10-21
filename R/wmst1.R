#############################
# wmst1 (one-arm) -- hidden
#############################



wmst1 <- function(time, status, tau0 = NULL, tau1 = NULL, alpha = 0.05) {

  #-- time
  #-- status
  #-- tau0 -- lower truncation time
  #-- tau1 -- upper truncation time
  #-- alpha -- gives (1-alpha) confidence interval

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
  if(alpha <=0 | alpha >= 1){
    stop(paste("The choice of alpha must be greater than 0 and less than 1."))
  }

  # compute window mean

  ft = survfit(Surv(time, status) ~ 1)
  df <- with(ft, data.frame(time, n.risk, n.event, surv))
  if (sum(which(df$n.event == 0)) != 0) {
    df <- df[-which(df$n.event == 0),]
  }

  df$term <- ifelse((df$n.risk - df$n.event) == 0, 0, df$n.event / (df$n.risk * (df$n.risk - df$n.event)))
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

  # create empy covariance matrix of correct dimensions

  if(tau0 < tau1){

    if (max(df$time) > tau0){
      t0 <- apply(X = outer(df$time, tau0, '>'),
                  MARGIN = 2, FUN = function(x) min(which(x)))
    } else {
      t0 <- which.max(df$time[which(df$time <= tau0)]) + 1
    }

    if (min(df$time) < tau1){
      t1 <- apply(X = outer(df$time, tau1, '<'), MARGIN = 2,
                  FUN = function(x) max(which(x))) + 1
    } else {
      t1 <- 1
    }

    cov.mat = matrix(0, nrow = length(t0:t1),length(t0:t1))

    # populate control group covariance matrix

    surv <- c(0,df$surv)
    term <- c(0, df$term)

    for (i in t0:t1) {
      for (j in t0:t1) {
        t <- min(i, j)
        S.ti <- surv[i]
        S.tj <- surv[j]
        risk.sum <- sum(term[1:t])
        cov.mat[i - t0 + 1, j - t0 + 1] <- S.ti * S.tj * risk.sum
      }
    }

    if(tail(time.diff, n = 1) == 0){
      time.diff <- head(time.diff, -1)
    }

    # compute variance

    time.diff <- as.matrix(time.diff)
    wmst.var <- t(time.diff) %*% cov.mat %*% time.diff
    wmst.se <- sqrt(wmst.var)

    # output

    out = matrix(0, 2, 4)
    out[1, ] = c(wmst, wmst.se, wmst - qnorm(1 - alpha/2) * wmst.se,
                 wmst + qnorm(1 - alpha/2) * wmst.se)
    out[2, ] = c((tau1 - tau0) - out[1, 1], wmst.se, (tau1 - tau0) - out[1, 4],
                 (tau1 - tau0) - out[1, 3])
    if (tau0 > 0) {
      rownames(out) = c("WMST", "WMTL")
    } else if (tau0 == 0) {
      rownames(out) = c("RMST", "RMTL")
    }
    colnames(out) = c("Est.", "se", paste("lower .",
                                          round((1 - alpha) * 100, digits = 0), sep = ""),
                      paste("upper .", round((1 - alpha) * 100, digits = 0),
                            sep = ""))
    Z = list()
    Z$note = NOTE
    Z$result = out
    Z$wmst = out[1, ]
    Z$wmtl = out[2, ]
    Z$tau0 = tau0
    Z$tau1 = tau1
    Z$wmst.var = wmst.var
    Z$fit = ft

  } else if (tau0 == tau1) {
    wmst.var <- (df$surv[idx]^2) * sum(df$term[1:idx])
    wmst.se <- sqrt(wmst.var)

    out = matrix(0, 1, 4)
    out[1,] = c(wmst,
                wmst.se,
                wmst - qnorm(1 - alpha / 2) * wmst.se,
                wmst + qnorm(1 - alpha / 2) * wmst.se)
    rownames(out) = c("Milestone")
    colnames(out) = c("Est.",
                      "se",
                      paste("lower .",
                            round((1 - alpha) * 100, digits = 0), sep = ""),
                      paste("upper .", round((1 - alpha) * 100, digits = 0),
                            sep = ""))
    Z = list()
    Z$note = NOTE
    Z$result = out
    Z$tau0 = tau0
    Z$tau1 = tau1
    Z$wmst.var = wmst.var
    Z$fit = ft

  }

  class(Z) = "wmst1"
  return(Z)


}

