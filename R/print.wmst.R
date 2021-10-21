#' @name print.wmst
#' @aliases print.wmst
#' @title print.wmst
#' @description S3 method for class 'wmst2'
#' @param x Object to be printed.
#' @param digits Integer indicating the number of decimal places.
#' @param ... Further arguments ignored in this function.
#' @return returns summary output for class 'wmst'
#' @export
######################################
# print.wmst (hidden)
######################################
print.wmst=function(x, digits=5, ...){

  cat("\n")

  cat(x$note,"\n\n")


    cat ("Window Mean Survival Time (WMST) by arm \n")

    prmatrix(round(x$wmst , digits=digits))

    cat("\n\n")

    cat ("Window Mean Time Lost (WMTL) by arm \n")

    prmatrix(round(x$wmtl, digits=digits))

    cat("\n\n")

    cat ("Between-group contrast \n")

    prmatrix(round(x$result, digits=digits))




  invisible(x)
}
NULL
