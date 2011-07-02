\name{regress}
\alias{regress}
\title{<title>}
\usage{regress(Y, X, family=binomial(), method="glm", formula=Y ~ ., ...)}
\description{<description>}
\details{<details>}
\value{<return>}
\author{Sam Lendle}
\arguments{\item{Y}{dependent variable for regression}
\item{X}{predictor variables for regression}
\item{family}{error distribution. Can be \code{gaussian} or \code{binomial}}
\item{method}{\code{"glm"} for glm or \code{"SL"} for SuperLearner}
\item{formula}{regression formula for glm. Dependent variable should be called \code{Y}.  Defaults to \code{Y~.} if not specified.}
\item{...}{additional options to be passed to SuperLearner}}
