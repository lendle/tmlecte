\name{tmle.cte}
\alias{tmle.cte}
\title{TMLE for statistical parameter \emph{E[E(Y|A=1,B)-E(Y|A=0,B)|A=a]}}
\usage{
  tmle.cte(Y, A, B, a = 0, Delta = NULL, Q = NULL,
    Qform = NULL,
    Q.SL.library = c("SL.glm", "SL.step", "SL.glm.interaction"),
    Qbounds = NULL, g1W = NULL, gform = NULL,
    gbound = 0.001, pDelta1 = NULL, g.Deltaform,
    g.SL.library = c("SL.glm", "SL.step", "SL.glm.interaction"),
    family = gaussian(),
    fluctuation = c("logistic", "linear"), alpha = 0.995,
    id = 1:length(Y), target = TRUE,
    tol = sqrt(.Machine$double.eps), maxiter = 100,
    verbose = FALSE)
}
\arguments{
  \item{Y}{continuous or binary outcome variable}

  \item{A}{binary treatment indicator, \code{1} -
  treatment, \code{0} - control}

  \item{B}{vector, matrix, or dataframe containing
  covariates to control for. Baseline covariates for the
  ATT, baseline covariates and mediator for NDEU}

  \item{a}{the treatment for which you want to estimate the
  conditional effect, \code{1} or \code{0}. \code{1}
  corresponds to the ATT, \code{0} corresponds to the ATN
  or NDE.}

  \item{Delta}{indicator of missing outcome or treatment
  assignment. \code{1} - observed, \code{0} - missing}

  \item{Q}{optional \emph{nx2} matrix of initial values for
  Q portion of the likelihood, \emph{(E(Y|A=0,B),
  E(Y|A=1,B))}}

  \item{Qform}{optional regression formula for estimation
  of \emph{E(Y|A,B)}, suitable for call to \code{glm}}

  \item{Q.SL.library}{optional vector of prediction
  algorithms to use for \code{SuperLearner} estimation of
  initial Q}

  \item{Qbounds}{vector of upper and lower bounds on Y and
  predicted values for initial \code{Q}. Defaults to the
  range of \code{Y}, widened by 10\% of the min and max
  values.}

  \item{g1W}{optional vector of conditional treatment
  assingment probabilities, \emph{P(A=1|B)}}

  \item{gform}{optional regression formula of the form A~W,
  if specified this overrides the call to
  \code{SuperLearner}}

  \item{gbound}{value between (0,1) for truncation of
  predicted probabilities. See \code{Details} section for
  more information}

  \item{pDelta1}{optional matrix of conditional
  probabilities for missingness mechanism with dimension
  \emph{nx2} \emph{P(Delta=1|A=0,B)},
  \emph{P(Delta=1|A=1,B)}}

  \item{g.Deltaform}{optional regression formula of the
  form \code{Delta~A+B}, if specified this overrides the
  call to \code{SuperLearner}}

  \item{g.SL.library}{optional vector of prediction
  algorithms to use for \code{SuperLearner} estimation of
  \code{g1W} or \code{pDelta1}}

  \item{family}{family specification for working regression
  models, generally 'gaussian' for continuous outcomes
  (default), 'binomial' for binary outcomes}

  \item{fluctuation}{'logistic' (default), or 'linear'}

  \item{alpha}{used to keep predicted initial values
  bounded away from (0,1) for logistic fluctuation}

  \item{id}{optional subject identifier}

  \item{target}{\code{TRUE} (default) to preform bias
  reduction (targeting) step. \code{FALSE} to return simple
  G computation based substitution estimate}

  \item{tol}{convergence criterion, defaults to square root
  of machine epsilon}

  \item{maxiter}{maximum number of iterations in the
  targeting step. Default is 100.}

  \item{verbose}{status messages printed if set to
  \code{TRUE} (default=\code{FALSE})}
}
\value{
  An object of class \code{cte} which is a list with the
  following components: \describe{
  \item{\code{psi}}{Parameter estimate}
  \item{\code{var.psi}}{Influence curve based estimate of
  the variance of \code{psi}} \item{\code{CI}}{95\%
  confidence interval for \code{psi}}
  \item{\code{pvalue}}{p-value from a two-sided test of
  \code{H0: psi_0=0}} \item{\code{IC}}{A vector with the
  value of the estimated influence curve for each
  observation} \item{\code{estimand}}{Character string
  describing the estimand, depends on the value of
  \code{a}} \item{\code{Qinit}}{Initial estimate of Q from
  the \code{estimateQ} function} \item{\code{Qstar}}{Matrix
  of final estimates of Q after targeting}
  \item{\code{ginit}}{Initital estimate of g from the
  \code{estimateg} function} \item{\code{g1Wstar}}{Vector
  of final estiates of g after targeting}
  \item{\code{g.Delta}}{Estimate of P(Delta=1|A,W) from the
  \code{estimateg} function} } If convergence fails, all
  items are either \code{NA} or \code{NULL}. If
  \code{target==FALSE}, some items are \code{NA} or
  \code{NULL}.
}
\description{
  Estimates the statistical parameter
  \emph{E[E(Y|A=1,B)-E(Y|A=0,B)|A=a]} using Targetted
  Maximum Likelihood.
}
\details{
  If \code{target} is \code{FALSE}, the targetting step is
  not done, and a G-computation type estimate based only on
  the inital estimate of Q is returned.  In this case, the
  influence curve, variance estimate, confidence interval
  and p-value are not calculated.

  Q.SL.library defaults to ('SL.glm', 'SL.step',
  'SL.glm.interaction') g.SL.library Defaults to ('SL.glm',
  'SL.step', 'SL.glm.interaction') This choice is simply
  because these algorithms are included in the base R
  installation. See SuperLearner help files for further
  information.
}
\author{
  Sam Lendle \email{lendle@stat.berkeley.edu}
}
\seealso{
  \code{\link{tmle.att}}, \code{\link{tmle.ndeu}},
  \code{\link{print.cte}}
}

