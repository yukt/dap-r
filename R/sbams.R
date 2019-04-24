#' Read a sbams file
#'
#' @param file file path to the sbams file
#' @param normalize TRUE by default. Data should be normalized before analysis.
#' @return a sbams object
#' @useDynLib dap, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @export
read.sbams <- function(file, normalize = TRUE) {
  .Call(`_dap_read_sbams`, PACKAGE = 'dap', file, normalize)
}

#' Integrative Genetic Association Analysis using Deterministic Approximation of Posteriors
#' 
#' The \code{dap.sbams} function is a simpler wrapper for \code{dap} function on sbams objects.
#'
#' @param sbams a sbams object
#' @return a dap object including model summary, SNP summary and signal cluster summary
#' @useDynLib dap, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @export
dap.sbams <- function(sbams)
{
  .Call(`_dap_dap`, PACKAGE = 'dap', list(data=sbams$file, output="output.dap"))
}

#' summary sbams object
#' 
#' \code{summary.sbams} helps to generate the summary statistics from a sbams object
#' to facilitate downstream analysis.
#'
#' @param sbams a sbams object
#' @return It will show the sample size [n] and total sum of squares for the quantitative trait [syy],
#' and write files [LD.dat] and [est.dat] under the working directory.
#' @useDynLib dap, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @export
summary.sbams <- function(sbams) {
  stopifnot(class(sbams)=="sbams")
  cat("============== SBAMS DATA SUMMARY ==============\n")
  .Call(`_dap_dap`, PACKAGE = 'dap', list(data=sbams$file, summary=TRUE))
  cat("Two summary files are generated:\n")
  cat("  [LD.dat]  - LD matrix between SNPs.\n")
  cat("  [est.dat] - effect estimates from single SNP testing.\n\n")
}


#' Integrative Genetic Association Analysis using Deterministic Approximation of Posteriors
#' 
#' The \code{dap} function accepts two different sets of inputs.
#' \enumerate{
#' \item sbams (either file path or object)
#' \item summary statistics (ld, est, n, syy)
#' }
#' @usage 
#' dap(file="/path/to/sbams/file") 
#' dap(sbams=a.sbams.object)
#' dap(ld="/path/to/ld/file", est="/path/to/est/file", n, syy)
#'
#' @param file  file path to a sbams file
#' @param sbams a sbams object
#' @param ld    file path to a ld file
#' @param est   file path to a est file
#' @param n     an integer indicating the sample size
#' @param syy   a float number indicating the total sum of squares for the quantitative trait
#' @return a dap object including model summary, SNP summary and signal cluster summary
#' 
#' 
#' @details The summary statistics can be obtained by \code{summary(a.sbams.object)}.
#' @useDynLib dap, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @export
dap <- function(file=NULL, sbams=NULL, ld=NULL, est=NULL, n=NA, syy=NA)
{
  if(!is.null(sbams))
  {
    stopifnot(class(sbams)=="sbams")
    return(.Call(`_dap_dap`, PACKAGE = 'dap', list(data=sbams$file, output="output.dap")))
  }
  if(!is.null(file))
  {
    return(.Call(`_dap_dap`, PACKAGE = 'dap', list(data=file, output="output.dap")))
  }
  if(!is.null(ld) & !is.null(est) & !is.na(n) & !is.na(syy))
  {
    return(.Call(`_dap_dap`, PACKAGE = 'dap', list(ld=ld,est=est,n=n,syy=syy, output="output.dap")))
  }
}