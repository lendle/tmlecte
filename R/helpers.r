##' .setColnames function
##'
##' assign names to every unnamed column of x
##' @title .setColnames function
##' @param x.colnames current column names
##' @param x.ncols current number of columns
##' @param firstChar prefix for internally assigned name
##' @return column names
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

##' .bound funciton 
##' 
##' set outliers to min/max allowable values
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

