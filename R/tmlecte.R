##' Prints results from tmle.cte.
##'
##' 
##' @title \code{print.cte} function
##' @param x an ojbect of class \code{cte}
##' @param ... ignored
##' @return Nothing
##' @author Sam Lendle
print.cte <- function(x, ...) {
  cat(x$estimand, ": ", x$psi,
      "\nEstimated Variance: ", x$var.psi,
      "\n95% Confidence interval: (", x$CI[1], ", ", x$CI[2],")",
      "\nP-value for H0: psi=0 vs. H1: psi!=0:", format.pval(x$pvalue, eps=0.001), "\n", sep="")
}

##' Estimates the ``conditional treatment effect'' E(E(Y|A=1,W)-E(Y|A=0,W)|A=a) using Targetted Maximum Likelihood.
##'
##' <details>
##' @title tmle.cte function
##' @param A binary treatment indicator where 1 is treated and 0 is
##' not treated 
##' @param B a matrix or dataframe of covariates. 
##' @param Y outcome variable. Currently should be in [0, 1], but can
##' be continuous.
##' @param a the treatment for which you want to estimate the
##' conditional effect, 1 or 0. 1 corresponds to the ATT, 0
##' corresponds to the ATN or NDE.
##' @param Q.method method to estimate \code{Q(A,B)=E(Y|A,B)}. Should be \code{"glm"} for glm or \code{"SL"} for SuperLearner.  If SuperLearner is not installed, glm is automatically used.
##' @param Q.formula the formula for glm if glm is used to estimate \code{Q(A,B)}. The response variable should be calle d \code{Y} and the predictors should include \code{A} and variables in \code{B}.  If missing, defaults to $Y~.$ (main terms including \code{A} all variables in \code{B}).
##' @param Q.SL.library SuperLearner library for estimating \code{Q(A,B)} if SuperLearner is used.
##' @param g.method method to estimate \code{g(A=a|B)=P(A=a|B)}. Should be \code{"glm"} for glm or \code{"SL"} for SuperLearner.  If SuperLearner is not installed, glm is automatically used.
##' @param g.formula the formula for glm if glm is used to estimate \code{g(A=a|B)}. The response variable should be calle d \code{A} and the predictors should include variables in \code{B}.  If missing, defaults to $A~.$ (main terms including all variables in \code{B}).
##' @param g.SL.library SuperLearner library for estimating \code{g(A=a|B)} if SuperLearner is used.
##' @param family \code{gaussian} for a continuous outcome \code{Y}, \code{binomial} for a
##' biniomial outcome
##' @param tol convergence criterion, set to something small
##' @param maxiter maximum number of targetting steps. 100 is good. If it takes more than that it's probably not going to converge.
##' @param target TRUE to preform bias reduction (targeting)
##' step. FALSE to return simple G computation based substitution
##' estimate
##' @param verbose TRUE for details for each iteration
##' @param Qbound a numeric vector containing minimum and maximum
##' bounds for Q. Should be between 0 and 1 for now.  Set to something close to 0 and 1. 
##' @param gbound a numeric vector containing minimum and maximum
##' bounds for g. Should be between 0 and 1. Set to something close to 0 and 1.
##' @param ... aditional parameters to be passed to (both)
##' SuperLearner calls 
##' @return An object of class \code{cte}
##' <details on returned object>
##' @author Sam Lendle \email{lendle@@stat.berkeley.edu}
##' @export
tmle.cte <- function(A, B, Y, a=0, Q.method="glm", Q.formula=NULL, Q.SL.library=NULL, g.method="glm", g.formula=NULL, g.SL.library=c("SL.glm", "SL.step", "SL.knn"), family=gaussian(), tol=1e-10, maxiter=100, target=TRUE, verbose=FALSE, Qbound=c(1e-10, 1-1e-10), gbound=c(1e-10, 1-1e-10), ...) {

  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  if (!(family$family %in% c("gaussian", "binomial"))) {
    stop("Currently only gaussian and binomial families are supported")
  }
  if (is.null(Q.SL.library)) {
    if (family$family=="gaussian") {
      Q.SL.library <- c("SL.glm", "SL.step")
    } else {
      Q.SL.library <- c("SL.glm", "SL.step", "SL.knn")
    }
  }

    
  Aa <- as.numeric(A==a)

  #Q.init.fit <- SuperLearner(Y, data.frame(A, B), SL.library=Q.SL.library, family=family, ...)
  Q.init.fit <- regress(Y, data.frame(A=A, B), family=family,
                        method=Q.method,
                        formula=Q.formula,
                        SL.library=Q.SL.library,
                        ...)
  #Q.init.fit <- glm(Qform, gaussian, dat, subset=(Delta==1))
  #Change to SL fit if you want, but remember that gform only has one variable so some
  #algorithms won't work
  Q.A1 <- predict(Q.init.fit, newdata=data.frame(A=1, B), X=data.frame(A=A, B), Y=Y)
  Q.A0 <- predict(Q.init.fit, newdata=data.frame(A=0, B), X=data.frame(A=A, B), Y=Y)

  
  
  #g.init.fit <- SuperLearner(A, B, SL.library=g.SL.library, family=binomial,...)
  g.init.fit <- regress(A, B, family=binomial,
                        method=g.method,
                        formula=g.formula,
                        SL.library=g.SL.library,
                        ...)
  g.A1 <- .bound(predict(g.init.fit), gbound)
  g.A0 <- 1-g.A1
  g.Aa <- a*g.A1 + (1-a)*g.A0

  
  fail=FALSE
  if (verbose) {
    cat("Variables in B:", names(B), "\n")
    if (!target) {
      cat("No targetting\n")
    } else {
      cat("Tolerance:", tol, "\nMaximum iterations:", maxiter, "\n")
    }
  }

  if (target) {
    iter=0
    done=FALSE
    prev.crit <- -1

    while(!done) {
      iter=iter+1

      #bound Q away from 0 and 1 so logit(Q) can be calculated
      Q.A1 <- .bound(Q.A1, Qbound)
      Q.A0 <- .bound(Q.A0, Qbound)

      #Calculate target parameter at current step
      psi <- mean(Q.A1[A==a] - Q.A0[A==a])
 
      #Calculate expected outcome under treatment recieved
      Q.AA <- A*Q.A1 + (1-A)*Q.A0

      #bound g away from 0 and 1 so logit(g) can be calculated, and clever covars
      g.Aa <- .bound(g.Aa, gbound)
      if (a==1) {
        g.A1 <- g.Aa
        g.A0 <- 1-g.A1
      } else {
        g.A0 <- g.Aa
        g.A1 <- 1-g.A0
      }
      #g.A1 <- a*g.Aa + (1-a)*(1-g.Aa)
      #g.A0 <- 1-g.A1
      
      
      #Calculate covariates for fluctuation
      #H1 <- (A*g.A0/(1-g.A0)-(1-A))
      #H1 <- A*g.Aa/g.A1 - (1-A)*g.Aa/g.A0
      
      #H1.A1 <- Delta/g.D1
      #H1.A0 <- (Delta/g.D1)*(-g.A1/(1-g.A1))
      H1.A1 <- g.Aa/g.A1
      H1.A0 <- -g.Aa/g.A0
      H1 <- A*H1.A1 + (1-A)*H1.A0
      
      H2 <- (Q.A1-Q.A0 - psi)

      #Estimate epsilon for logistic fluctuation of Q and g
      Q.up <- suppressWarnings(glm(Y~-1+H1, offset=qlogis(Q.AA), family=binomial))
      #Q.up <- glm(Y~-1+H1, dat, offset=Q.AA, family=gaussian, subset=(Delta==1))
      g.up <- glm(Aa~-1+H2, offset=qlogis(g.Aa), family=binomial)
      eps <- c(coef(Q.up), coef(g.up))
      if (is.na(eps[2])) {eps[2] <- 0}

      #Fluctuate predicted outcome and treatment, Q and g
      Q.A1 <- plogis(qlogis(Q.A1)+eps[1]*H1.A1)
      #Q.A1 <- Q.A1+eps[1]*H1.A1
      Q.A0 <- plogis(qlogis(Q.A0)+eps[1]*H1.A0)
      #Q.A0 <- Q.A0+eps[1]*H1.A0
      g.Aa <- plogis(qlogis(g.Aa)+eps[2]*H2)
      
      crit <- max(abs(na.omit(eps)))
      if (crit <= tol || iter>=maxiter) {
        done=TRUE
      }

      #If not done and critical value has not improved (probably because of bounding on Q)
      #then convergence failed
      #Note that if Q is not bounded, logit(Q)=+/-Inf so convergence fails either way...
      if(done==FALSE && crit==prev.crit) {
        done=TRUE
        fail=TRUE
        warning("Convergence criterion failed to improve between iterations")
      }
      prev.crit <- crit
      
      if (verbose) cat("iter:", iter, "eps:", eps, "crit:", crit, "psi:", psi, "\n")
    }
  }

  if (fail) {
    res <- list(psi=NA,
                var.psi=NA,
                CI=NA,
                pvalue=NA
                )
  }

  if (!fail) {
    #Calculate final parameter of interest after convergence
    psi <-  mean(Q.A1[A==a] - Q.A0[A==a])
    Q.AA <- A*Q.A1 + (1-A)*Q.A0
    p.Aa <- mean(Aa)
    IC <- ((A*g.Aa/g.A1 - (1-A)*g.Aa/g.A0)*(Y-Q.AA)+Aa*(Q.A1-Q.A0 - psi))/p.Aa
    if (verbose) {
      cat("Mean of the influence curve:", mean(IC), "\n")
    }
    var.psi <- var(IC)/length(Y)
    CI <- psi + c(-1.96, 1.96)*sqrt(var.psi)
    pvalue <- 2*(1-pnorm(abs(psi/sqrt(var.psi))))
    estimand <- paste("E(E(Y|A=1,B)-E(Y|A=0,B)|A=", a, ")", sep="")
    res <- list(psi=psi, var.psi=var.psi, CI=CI, pvalue=pvalue, IC=IC, estimand=estimand)
  }
  class(res) <- c("cte")
  res
}

