################################
# tau.pairs (one-arm) -- hidden
################################

tau.pairs <- function(time = NULL, status = NULL, type = NULL, round = FALSE, combination = FALSE, tau1_max = NULL){
  
  
  if (combination) {
    if (type == "quarter" & round) {
      tau0 = round(seq(from = 0, to = 1, by = 0.25)[1:4] * tau1_max)
      tau1 = c(round(seq(from = 0, to = 1, by = 0.25)[2:4] * tau1_max), tau1_max)
    } else if (type == "quarter" & !round) {
      tau0 = seq(from = 0, to = 1, by = 0.25)[1:4] * tau1_max
      tau1 = c(seq(from = 0, to = 1, by = 0.25)[2:4] * tau1_max, tau1_max)
    } else if (type == "quantile" & round) {
      tau0 = round(seq(from = 0, to = 1, by = 0.1)[1:10] * tau1_max)
      tau1 = c(round(seq(from = 0, to = 1, by = 0.1)[2:10] * tau1_max), tau1_max)
    } else if (type == "quantile" & !round) {
      tau0 = seq(from = 0, to = 1, by = 0.1)[1:10] * tau1_max
      tau1 = c(seq(from = 0, to = 1, by = 0.1)[2:10] * tau1_max, tau1_max)
    } else if (type == "quartile" & round){
      tau0 = c(0,round(quantile(time[status == 1])[2:4]))
      tau1 = c(round(quantile(time[status == 1])[2:4]),tau1_max)
    } else if (type == "quartile" & !round){
      tau0 = c(0,quantile(time[status == 1])[2:4])
      tau1 = c(quantile(time[status == 1])[2:4],tau1_max)
    } else if (type == "survquart" & round){
      fit <- survfit(Surv(time,status)~1)
      quart <- (1 - tail(fit$surv,1))/4
      survquart <- c(fit$time[which.min(abs(fit$surv-(1 - quart)))],fit$time[which.min(abs(fit$surv-(1 - 2*quart)))],fit$time[which.min(abs(fit$surv-(1 - 3*quart)))])
      tau0 = c(0,round(survquart))
      tau1 = c(round(survquart),tau1_max)
    } else if (type == "survquart" & !round){
      fit <- survfit(Surv(time,status)~1)
      quart <- (1 - tail(fit$surv,1))/4
      survquart <- c(fit$time[which.min(abs(fit$surv-(1 - quart)))],fit$time[which.min(abs(fit$surv-(1 - 2*quart)))],fit$time[which.min(abs(fit$surv-(1 - 3*quart)))])
      tau0 = c(0,survquart)
      tau1 = c(survquart,tau1_max)
    }
    if (any(tau0 == tau1)) {
      tau0[which(tau0 == tau1) + 1] <- tau0[which(tau0 == tau1) + 1] + 1
      tau1[which(tau0 == tau1)] <- tau1[which(tau0 == tau1)] + 1
    }
    
    tau.grid = expand.grid(tau0, tau1)
    names(tau.grid) = c("tau0", "tau1")
    tau = tau.grid[(tau.grid$tau0 < tau.grid$tau1), ]
    
  } else if (!combination) {
    if (type == "quarter" & round) {
      tau0 = round(seq(from = 0, to = 1, by = 0.25)[1:4] * tau1_max)
      tau1 = c(round(seq(from = 0, to = 1, by = 0.25)[2:4] * tau1_max), tau1_max)
    } else if (type == "quarter" & !round) {
      tau0 = seq(from = 0, to = 1, by = 0.25)[1:4] * tau1_max
      tau1 = c(seq(from = 0, to = 1, by = 0.25)[2:4] * tau1_max, tau1_max)
    } else if (type == "quantile" & round) {
      tau0 = round(seq(from = 0, to = 1, by = 0.1)[1:10] * tau1_max)
      tau1 = c(round(seq(from = 0, to = 1, by = 0.1)[2:10] * tau1_max), tau1_max)
    } else if (type == "quantile" & !round) {
      tau0 = seq(from = 0, to = 1, by = 0.1)[1:10] * tau1_max
      tau1 = c(seq(from = 0, to = 1, by = 0.1)[2:10] * tau1_max, tau1_max)
    } else if (type == "quartile" & round){
      tau0 = c(0,round(quantile(time[status == 1])[2:4]))
      tau1 = c(round(quantile(time[status == 1])[2:4]),tau1_max)
    } else if (type == "quartile" & !round){
      tau0 = c(0,quantile(time[status == 1])[2:4])
      tau1 = c(quantile(time[status == 1])[2:4],tau1_max)
    } else if (type == "survquart" & round){
      fit <- survfit(Surv(time,status)~1)
      quart <- (1 - tail(fit$surv,1))/4
      survquart <- c(fit$time[which.min(abs(fit$surv-(1 - quart)))],fit$time[which.min(abs(fit$surv-(1 - 2*quart)))],fit$time[which.min(abs(fit$surv-(1 - 3*quart)))])
      tau0 = c(0,round(survquart))
      tau1 = c(round(survquart),tau1_max)
    } else if (type == "survquart" & !round){
      fit <- survfit(Surv(time,status)~1)
      quart <- (1 - tail(fit$surv,1))/4
      survquart <- c(fit$time[which.min(abs(fit$surv-(1 - quart)))],fit$time[which.min(abs(fit$surv-(1 - 2*quart)))],fit$time[which.min(abs(fit$surv-(1 - 3*quart)))])
      tau0 = c(0,survquart)
      tau1 = c(survquart,tau1_max)
    }
    if (any(tau0 == tau1)) {
      tau0[which(tau0 == tau1) + 1] <- tau0[which(tau0 == tau1) + 1] + 1
      tau1[which(tau0 == tau1)] <- tau1[which(tau0 == tau1)] + 1
    }
    tau = data.frame("tau0" = tau0, "tau1" = tau1)
  }
  
  return(tau)
  
}