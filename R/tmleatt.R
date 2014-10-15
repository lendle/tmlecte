##' TMLE for the average treatment effect among the treated (ATT)
##'
##' This function estimates E[E(Y|A=1, W) - E(Y|A=0, W) | A=1].
##' Under particular causal assumptions, this is equal to the ATT.
##' This is a wrapper function to \code{tmle.cte}.
##' @title tmle.att function
##' @param A treatment
##' @param W baseline covariates
##' @param Y outcome
##' @param ... other arguments to be passed to \code{tmle.cte}
##' @return an object of class "cte"
##' @author Sam Lendle
##' @export
tmle.att <- function(A, W, Y, ...) {
  t <- tmle.cte(A, W, Y, a=1, ...)
  t$estimand <- "ATT"
  return(t)
}

