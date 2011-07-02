#---------- function .setColnames ---------------
# assign names to every unnamed column of x
# arguments
# 	x.colnames - current column names
#	x.ncols - current number of columns
# 	firstChar - prefix for internally assigned name
# return the names
#-----------------------------------------
##' <description>
##'
##' <details>
##' @title <title>
##' @param x.colnames <description>
##' @param x.ncols <description>
##' @param firstChar <description>
##' @return <return>
##' @author Susan Gruber
##' @name setColnames
##' @aliases .setColnames
.setColnames <- function(x.colnames, x.ncols, firstChar){
	if(is.null(x.colnames)) {
		if(x.ncols > 1){
			x.colnames <- paste(firstChar,1:x.ncols, sep="")
		} else {
			x.colnames <- firstChar
		}
	} else {
		invalid.name <- nchar(x.colnames) == 0
		if(any(invalid.name)){
			x.colnames[invalid.name] <- paste(".internal",firstChar, which(invalid.name), sep="")
		}
	}
	return(x.colnames)
}

##' set outliers to min/max allowable values
##'
##' @title .bound function
##' @param x object of numerical data, usually a matrix or vector
##' @param bounds object of numeric data where the min is the minimum
##' bound and the max is the maximum bound
##' @return x with values less than \code{min(bounds)} set to \code{min(bounds)}, likewise for values greater than \code{max(bounds)}
##' @author Susan Gruber
##' @name bound
##' @aliases .bound
.bound <- function(x, bounds){
	x[x>max(bounds)] <- max(bounds)
	x[x<min(bounds)] <- min(bounds)
	return(x)
}

##' <description>
##'
##' <details>
##' @title gendata function
##' @param n asdf
##' @param include.D aasdf
##' @param A.coef asdf
##' @param pDelta asdf
##' @return a data frame
##' @author Sam Lendle
##' @export
gendata <- function(n, include.D=FALSE, A.coef=0, pDelta=NULL) {
  if (is.null(pDelta)) pDelta <- function(...){1}
  #If include.D is false, D will be set to one for everyone
  if (include.D) {
    D <- rbinom(n, 1, 0.7) + 1
  } else {
    D <- rep(1, n)
  }
  W1 <- rbinom(n, 1, 0.35)
  W2 <- rbinom(n, 1, 0.6-0.3*(D-1))
  W3 <- rnorm(n, 1, D)
  W4 <- rnorm(n, D-1, 1)
  pa <- numeric(n)
  #pa[D==1] <- logistic(-1+1.7*W1+0.25*W3)[D==1]
  pa[D==1] <- plogis(-1+1.7*W1+0.25*W3+W4)[D==1]
  #pa[D==1] <- logistic(-3+1.7*W1+0.25*W3+4*(W4>=0))[D==1]
  pa[D==2] <- plogis(-1+W1-0.2*W2 + 0.4*W4)[D==2]
  A <- rbinom(n, 1, pa)
  Y <- rbinom(n, 1, plogis(-1+A.coef*A+0.5*W1-2*W2))
  Delta <- rbinom(n, 1, pDelta(A=A, W1=W1, W2=W2, W3=W3, W4=W4))
  data.frame(D, W1, W2, W3, W4, A, Y, Delta, pa)
}
