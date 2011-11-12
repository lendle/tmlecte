\name{tmle.cte}
\alias{tmle.cte}
\title{tmle.cte function}
\usage{
  tmle.cte(A, B, Y, a = 0, Delta = NULL, Q.method = "glm",
    Q.formula = NULL, Q.SL.library = NULL,
    g.method = "glm", g.formula = NULL,
    g.SL.library = c("SL.glm", "SL.step", "SL.knn"),
    gDelta.method = "glm", gDelta.formula = NULL,
    gDelta.SL.library = c("SL.glm", "SL.step", "SL.knn"),
    family = gaussian(), tol = 1e-10, maxiter = 100,
    target = TRUE, verbose = FALSE,
    Qbound = c(1e-10, 1 - 1e-10),
    gbound = c(1e-10, 1 - 1e-10), ...)
}
\arguments{
  \item{A}{binary treatment indicator where 1 is treated
  and 0 is not treated}

  \item{B}{a matrix or dataframe of covariates.}

  \item{Y}{outcome variable. Currently should be in [0, 1],
  but can be continuous.}

  \item{a}{the treatment for which you want to estimate the
  conditional effect, 1 or 0. 1 corresponds to the ATT, 0
  corresponds to the ATN or NDE.}

  \item{Delta}{binary indicator of missing outcome.  Should
  be 0 where \code{Y} is \code{NA} or should be treated as
  \code{NA}, and 1 otherwise.}

  \item{Q.method}{method to estimate
  \code{Q(A,B)=E(Y|A,B)}. Should be \code{"glm"} for glm or
  \code{"SL"} for SuperLearner.  If SuperLearner is not
  installed, glm is automatically used.}

  \item{Q.formula}{the formula for glm if glm is used to
  estimate \code{Q(A,B)}. The response variable should be
  calle d \code{Y} and the predictors should include
  \code{A} and variables in \code{B}.  If missing, defaults
  to \code{Y~.} (main terms including \code{A} all
  variables in \code{B}).}

  \item{Q.SL.library}{SuperLearner library for estimating
  \code{Q(A,B)} if SuperLearner is used.}

  \item{g.method}{method to estimate
  \code{g(A=a|B)=P(A=a|B)}. Should be \code{"glm"} for glm
  or \code{"SL"} for SuperLearner.  If SuperLearner is not
  installed, glm is automatically used.}

  \item{g.formula}{the formula for glm if glm is used to
  estimate \code{g(A=a|B)}. The response variable should be
  calle d \code{A} and the predictors should include
  variables in \code{B}.  If missing, defaults to $A~.$
  (main terms including all variables in \code{B}).}

  \item{g.SL.library}{SuperLearner library for estimating
  \code{g(A=a|B)} if SuperLearner is used.}

  \item{gDelta.method}{method to estimate
  \code{P(Delta=1|A, B)}. Should be \code{"glm"} for glm or
  \code{"SL"} for SuperLearner.  If SuperLearner is not
  installed, glm is automatically used.}

  \item{gDelta.formula}{the formula for glm if glm is used
  to estimate \code{P(Delta=1|A, B)}. The response variable
  should be called \code{Delta} and the predictors should
  include \code{A} and variables in \code{B}.  If missing,
  defaults to $Delta~.$ (main terms including \code{A} and
  all variables in \code{B}).}

  \item{gDelta.SL.library}{SuperLearner library for
  estimating \code{P(Delta=1|A,B)} if SuperLearner is
  used.}

  \item{family}{\code{gaussian} for a continuous outcome
  \code{Y}, \code{binomial} for a biniomial outcome}

  \item{tol}{convergence criterion, set to something small}

  \item{maxiter}{maximum number of targetting steps. 100 is
  good. If it takes more than that it's probably not going
  to converge.}

  \item{target}{TRUE to preform bias reduction (targeting)
  step. FALSE to return simple G computation based
  substitution estimate}

  \item{verbose}{TRUE for details for each iteration}

  \item{Qbound}{a numeric vector containing minimum and
  maximum bounds for Q. Should be between 0 and 1 for now.
  Set to something close to 0 and 1.}

  \item{gbound}{a numeric vector containing minimum and
  maximum bounds for g. Should be between 0 and 1. Set to
  something close to 0 and 1.}

  \item{...}{aditional parameters to be passed to (both)
  SuperLearner calls}
}
\value{
  An object of class \code{cte} <details on returned
  object>
}
\description{
  Estimates the ``conditional treatment effect''
  E(E(Y|A=1,W)-E(Y|A=0,W)|A=a) using Targetted Maximum
  Likelihood.
}
\details{
  <details>
}
\author{
  Sam Lendle \email{lendle@stat.berkeley.edu}
}

