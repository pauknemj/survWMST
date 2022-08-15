#' @name print.pow
#' @aliases print.pow
#' @title print.pow
#' @description S3 method for class 'pow'
#' @param x Object to be printed.
#' @param digits Integer indicating the number of decimal places.
#' @param ... Further arguments ignored in this function.
#' @return returns summary output for class 'pow'
#' @export
######################################
# print.pow (hidden)
######################################
print.pow=function(x, digits=5, ...){
  
  cat("\n")
  
  if(x$two.sided == TRUE){
    cat(past0("The following output is for an two-sided WMST test with window bounds \n(tau0,tau1) = (",x$tau0,",",x$tau1,")"," based on ",x$NOTE,": \n\n"))
  } else if(x$two.sided == FALSE){
    cat(paste0("The following output is for an upper-tailed WMST test with window bounds \n(tau0,tau1) = (",x$tau0,",",x$tau1,")"," based on ",x$NOTE,": \n\n"))
  }
  
  if(x$indicator == "power"){
    
    cat("The total required sample size: \n")
    print(x$n)
    
    cat("\n")
    
    cat(paste0("The actual asymptotic power estimate based on a sample size of ", x$n,": \n"))
    print(x$power)
    
  } else if(x$indicator == "n"){
    
    cat("The power of the test: \n")
    print(x$power)
    
  }
  
  

  invisible(x)
}
NULL
