##' <description>
##'
##' <details>
##' @title <title>
##' @param Y continuous or binary outcome variable
##' @param A binary treatment indicator, \code{1} - treatment, \code{0} - control
##' @param WZ vector, matrix, or dataframe containing baseline covariate(s), if applicable, and mediator(s).
##' @param RCT \code{TRUE} indicates an RCT, wich changes the label of the estimand in the printed output to "Natural direct effect".
##' @param ... other arguments to be passed to \code{\link{tmle.cte}}.
##' @return an object of class \code{cte} with \code{$estimatd} set to "Natural direct effect (among the untreated)".
##' @seealso \code{\link{tmle.cte}}
##' @author Sam Lendle \email{lendle@@stat.berkeley.edu}
##' @export
tmle.ndeu <- function(Y, A, WZ, RCT=FALSE, ...) {
  t <- tmle.cte(Y=Y, A=A, B=WZ, a=0, ...)
  t$estimand <- ifelse(RCT,
                       "Natural direct effect",
                       "Natural direct effect among the untreated")
  return(t)
}

