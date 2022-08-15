#' @name print.combo
#' @aliases print.combo
#' @title print.combo
#' @description S3 method for class 'combo'
#' @param x Object to be printed.
#' @param digits Integer indicating the number of decimal places.
#' @param ... Further arguments ignored in this function.
#' @return returns summary output for class 'combo'
#' @export
######################################
# print.combo (hidden)
######################################
print.combo=function(x, digits=5, ...){
  
  cat("\n")
  
  cat(x$NOTE,"\n\n")
  
  
  cat ("P-value from versatile WMST with method =", x$method,"and alternative =", x$alternative,": \n")
  
  print(round(x$p , digits=7))
  
  cat("\n\n")
  
  cat ("Difference in WMST for each window: \n")
  
  prmatrix(round(x$result, digits=digits))
  
  
  
  
  invisible(x)
}
NULL
