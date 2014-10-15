##' TMLE for the NDE in a randomized control trial or for the NDE among the untreated
##'
##' Estimates the statistical parameter E(E(Y|A=1,W,Z)-E(Y|A=0,W,Z)|A=0)
##' Under particular causal assuptions, this can be interpreted as the natural direct effect (NDE) 
##' among the untreated.
##' When A is indepdent of W as in an RCT, this is equal to the NDE. 
##' This is a wrapper function to \code{tmle.cte}.
##'  
##' @title tmle.nde function
##' @param A treatment
##' @param WZ baseline covariates and mediators
##' @param Y outcome
##' @param RCT TRUE indicates an RCT. This changes the label of the estimand in the output
##' from "Natural direct effect among the untreated" to "Natural direct effect".
##' @param ... other arguments to be passed to \code{tmle.cte}
##' @return an object of class "cte"
##' @author Sam Lendle
##' @seealso tmle.cte
##' @export
tmle.nde <- function(A, WZ, Y, RCT=FALSE, ...) {
  t <- tmle.cte(A, WZ, Y, a=0, ...)
  t$estimand <- ifelse(RCT,
                       "Natural direct effect",
                       "Natural direct effect among the untreated")
  return(t)
}

