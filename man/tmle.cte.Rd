\name{tmle.cte}
\alias{tmle.cte}
\title{tmle.cte function}
\usage{tmle.cte(A, B, Y, a=0, Q.method="glm", Q.formula, Q.SL.library,
    g.method="glm", g.formula, g.SL.library=c("SL.glm", "SL.step",
    "SL.knn"), family=gaussian(), tol=1e-10, maxiter=100, target=TRUE,
    verbose=FALSE, Qbound=c(1e-10, 1 - 1e-10), gbound=c(1e-10, 1 -
    1e-10), ...)}
\description{Estimates the ``conditional treatment effect'' E(E(Y|A=1,W)-E(Y|A=0,W)|A=a) using Targetted Maximum Likelihood.}
\details{<details>}
\value{adsf}
\author{Sam Lendle \email{lendle@stat.berkeley.edu}}
\arguments{\item{A}{binary treatment indicator where 1 is treated and 0 is
not treated}
\item{B}{a matrix or dataframe of covariates. In an RCT,
can include baseline covariates here to control for baseline
confounding}
\item{Y}{outcome variable. Currently should be in [0, 1], but can
be continuous.}
\item{a}{the treatment for which you want to estimate the
conditional effect, 1 or 0. 1 corresponds to the ATT, 0
corresponds to the ATN or NDE.}
\item{Q.SL.library}{SuperLearner library for estimating Q}
\item{g.SL.library}{SuperLearner library for estimating g}
\item{family}{gaussian for a continuous outcome, binomial for a
biniomial outcome}
\item{tol}{convergence criterion, set to something small}
\item{maxiter}{maximum number of targetting steps}
\item{target}{TRUE to preform bias reduction step. FALSE to return
simple G computation based estimate}
\item{verbose}{TRUE for details for each iteration}
\item{Qbound}{bound Q away from 0 and 1}
\item{gbound}{bound g away from 0 and 1}
\item{...}{aditional parameters to be passed to (both)
SuperLearner calls}}
