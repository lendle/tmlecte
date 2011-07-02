\name{bounds}
\title{.bound function}
\usage{.bound(x, bounds)}
\description{set outliers to min/max allowable values}
\value{x with values less than \code{min(bounds)} set to \code{min(bounds)}, likewise for values greater than \code{max(bounds)}}
\author{Susan Gruber}
\alias{.bounds}
\arguments{\item{x}{object of numerical data, usually a matrix or vector}
\item{bounds}{object of numeric data where the min is the minimum bound and the max is the maximum bound}}
