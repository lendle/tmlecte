##' Estimates the \strong{A}verage \strong{T}reatment effect among the \strong{T}reated (ATT) or \strong{U}ntreated (ATU), \emph{E[Y(1)-Y(O)|A=1]} or \emph{E[Y(1)-Y(0)|A=0]}, using targeted maximum likelihood.
##'
##' 
##' @title TMLE for the ATT or ATU
##' @param Y continuous or binary outcome variable
##' @param A binary treatment indicator, \code{1} - treatment, \code{0} - control
##' @param W vector, matrix, or dataframe containing baseline covariates to control for.
##' @param param The parameter being estimated, "ATT" (default) or "ATU".
##' @param ... other arguments passed to \code{\link{tmle.cte}}
##' @return An object of class \code{cte} with \code{$estimand} set to \code{param}.
##' @seealso \code{\link{tmle.cte}}
##' @author Sam Lendle \email{lendle@@stat.berkeley.edu}
tmle.att <- function(Y, A, W, param=c("ATT", "ATU"), ...) {
  param <- match.arg(param)
  if (param=="ATT") {
    a <- 1
  } else {
    a <- 0
  }
  t <- tmle.cte(Y=Y, A=A, B=W, a=a, ...)
  t$estimand <- param
  return(t)
}

