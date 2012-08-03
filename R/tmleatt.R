##' <description>
##'
##' <details>
##' @title <title>
##' @param Y outcome
##' @param A treatment
##' @param W baseline covariates
##' @param ... other arguments to be passed to tmle.cte
##' @return an object of class "cte"
##' @author Sam Lendle
##' @export
tmle.att <- function(Y, A, W, ...) {
  t <- tmle.cte(Y=Y, A=A, B=W, a=1, ...)
  t$estimand <- "ATT"
  return(t)
}

