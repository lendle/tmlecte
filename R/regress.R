##' <description>
##'
##' <details>
##' @title <title>
##' @param resp dependent variable for regression
##' @param X predictor variables for regression
##' @param family error distribution. Can be \code{gaussian} or \code{binomial}
##' @param method \code{"glm"} for glm or \code{"SL"} for SuperLearner
##' @param formula regression formula for glm. Defaults to \code{resp~.} if not specified.  The dependent variable on the lhs of the ~ does not matter, and is automatically changed to \code{resp}.
##' @param ... additional options to be passed to SuperLearner
##' @return <return>
##' @author Sam Lendle
##' @export
regress <- function(resp, X, family=binomial(), method="glm", formula=resp ~ ., ...) {
  SL.installed <- "SuperLearner" %in% installed.packages()[,1]
  SL.version <- NULL
  if (method=="glm" || !SL.installed) {
    if (method=="SL") {
      warning("SuperLearner is not installed, using main terms glm", call.=FALSE)
      method <- "glm"
    }
    if (is.null(formula)) {
      formula <- resp ~ .
    } else {
      formula[[2]] <- as.name("resp") # changes the lhs of the formula to resp
    }
    fit <- glm(formula, data=data.frame(resp, X), family=family)
  } else if (method=="SL") {
    require(SuperLearner)
    SL.version <- packageVersion("SuperLearner")$major
    if (SL.version==1) {
      warning("Your version of SuperLearner is out of date. You should consider updating to version 2 if you don't have a good reason not to...")
      fit <- SuperLearner(resp, data.frame(X), family=family, ...)
    }
    else {
      fit <- SuperLearner(resp, data.frame(X), family=family, ...)
    }
  }
  res <- list(fit=fit, method=method, SL.version=SL.version)
  class(res) <- "regress"
  return(res)
}

##' <description>
##'
##' <details>
##' @title <title>
##' @param object an object of class \code{regress}
##' @param newdata optional new data for prediction in the same format
##' as \code{X} used in the original \code{regress} call
##' @param X original predictors from the \code{regress} call. Needed
##' for some prediction algorithms when using \code{newdata}
##' @param Y original dependent variable from \code{regress}
##' call. Needed for some prediction algorithms when using
##' \code{newdata}
##' @param ... ignored
##' @return <return>
##' @author Sam Lendle
##' @method predict regress
##' @export
predict.regress <- function(object, newdata, X=NULL, Y=NULL, ...) {
  if (object$method=="glm") {
    if (missing(newdata)) return(predict(object$fit, type="response"))
    return(predict(object$fit, newdata=newdata, type="response"))
  }
  else if (object$method=="SL") {
    if (any(is.null(Y), is.null(X)) & !missing(newdata)) warning("Original data needs to be passed to predict.regress when using SuperLearner and newdata.  predict may fail depending on the SL.library otherwise...")
    if (object$SL.version==1) {
      if (missing(newdata)) {
        return(predict(object$fit))
      } else {
        return(predict(object$fit, newdata=data.frame(newdata), X=data.frame(X), Y=Y)$fit)
      }
    } else {
      if(missing(newdata)) return(predict(object$fit)$pred)
      return(predict(object$fit, newdata=data.frame(newdata), X=data.frame(X), Y=Y)$pred)
    }
  }
}
