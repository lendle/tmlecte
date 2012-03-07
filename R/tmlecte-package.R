##' Estimates treatment effects averaged over a distribution of covariates conditional on treatment
##'
##' \tabular{ll}{
##' Package: \tab tmlecte\cr
##' Type: \tab Package\cr
##' Version: \tab 0.2-1\cr
##' Date: \tab 2011-06-29\cr
##' License: \tab GPL (>= 2)\cr
##' LazyLoad: \tab yes\cr
##' }
##'
##' Estimates the Average Treatment effect among the Treated (ATT), or among
##' the not treated. The latter is the same statistical parameter as the Natural Direct
##' Effect (NDE) when there is no baseline confounding for either the treatment, A, or
##' the mediator(s)/intermediate variable(s), Z. When there is/are (a) baseline
##' confounder(s) (W), the "NDE among the not treated" is estimated. In the case of a
##' randomized trial, this is equivalent to the NDE, because P(W|A)=P(W).  In this case,
##' the NDE can be efficiently estimated were both Z and W have large dimensions. This
##' is a working package.  The methods in this package will likely be made available in
##' the tmle package on CRAN in the future.
##' @name tmlecte-package
##' @aliases tmlecte
##' @title Estimates the NDE and ATT with TMLE
##' @author Sam Lendle \email{lendle@@stat.berkeley.edu}
##' @references <ref>
##' @keywords package
{}
