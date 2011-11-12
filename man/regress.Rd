\name{regress}
\alias{regress}
\title{<title>}
\usage{
  regress(resp, X, family = binomial(), method = "glm",
    formula = resp ~ ., ...)
}
\arguments{
  \item{resp}{dependent variable for regression}

  \item{X}{predictor variables for regression}

  \item{family}{error distribution. Can be \code{gaussian}
  or \code{binomial}}

  \item{method}{\code{"glm"} for glm or \code{"SL"} for
  SuperLearner}

  \item{formula}{regression formula for glm. Defaults to
  \code{resp~.} if not specified.  The dependent variable
  on the lhs of the ~ does not matter, and is automatically
  changed to \code{resp}.}

  \item{...}{additional options to be passed to
  SuperLearner}
}
\value{
  <return>
}
\description{
  <description>
}
\details{
  <details>
}
\author{
  Sam Lendle
}

