#' summary
#'
#' @param object a dap object
#' @export
summary.dap = function(object, ...)
{
  stopifnot(class(object) == "dap")

  # Print title
  cat("============== Posterior Inference {DAP} ==============\n")
  cat(" Posterior expected model size:", round(object$size.mean,3), "( sd =", round(object$size.sd,3),")\n")
  cat(" LogNC =", round(object$logNC,5), "( Log10NC =", round(object$log10NC,3),")\n")
  cat(" Minimum PIP is estimated at", format(object$PIP.min,scientific = TRUE), "( N =", round(object$N),")\n")
  cat("=======================================================\n")
}
