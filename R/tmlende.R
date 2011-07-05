tmle.nde <- function(A, WZ, Y, RCT=FALSE, ...) {
  t <- tmle.cte(A, WZ, Y, a=0, ...)
  t$estimand <- ifelse(RCT,
                       "Natural direct effect",
                       "Natural direct effect among the untreated")
  return(t)
}

