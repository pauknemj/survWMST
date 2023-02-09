#' Plots two survival distributions.
#'
#' @import ggplot2
#' @param survdefC the survival distribution of the control group (will be
#' plotted as a solid line), as a list in the form output by `survdef`.
#' @param survdefT the survival distribution of the control group (will be
#'  plotted as a dashed line), as a list in the form output by `survdef`.
#' @param xupper the upper x axis limit for the plot.
#'
#' @export

################################
# survplot
################################


survplot <- function (survC = NA, survT = NA, xupper = NA) {

  x <- seq(from = 0, to = xupper, length.out = 100)
  data <- data.frame("x" = x, "survC" = survC$S(x), "survT" = survT$S(x))
  ggplot(data) +
    geom_line(aes(x = x, y = survC), linetype = "solid") +
    geom_line(aes(x = x, y = survT), linetype = "dashed") +
    ylim(0,1) + xlim(0,xupper) + labs(x = "time", y = "survival probability") +
    theme_bw()

}
