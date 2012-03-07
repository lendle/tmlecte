##' <description>
##'
##' <details>
##' @title <title>
##' @param A treatment
##' @param W baseline covariates
##' @param Y outcome
##' @param ... other arguments to be passed to tmle.cte
##' @return an object of class "cte"
##' @author Sam Lendle
##' @export
tmle.att <- function(A, W, Y, ...) {
  t <- tmle.cte(A, W, Y, a=1, ...)
  t$estimand <- "ATT"
  return(t)
}

