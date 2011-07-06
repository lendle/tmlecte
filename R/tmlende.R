##' <description>
##'
##' <details>
##' @title <title>
##' @param A treatment
##' @param WZ baseline covariates and mediators
##' @param Y outcome
##' @param RCT TRUE indicates and RCT. This only changes how the parameter is labeled in the output...
##' @param ... other arguments to be passed to tmle.cte
##' @return an object of class "cte"
##' @author Sam Lendle
##' @export
tmle.nde <- function(A, WZ, Y, RCT=FALSE, ...) {
  t <- tmle.cte(A, WZ, Y, a=0, ...)
  t$estimand <- ifelse(RCT,
                       "Natural direct effect",
                       "Natural direct effect among the untreated")
  return(t)
}

