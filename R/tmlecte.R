logit <- function(x) qlogis(x)
logistic <- function(x) plogis(x)

print.cte <- function(x, ...) {
  cat(estimand, ": ", x$psi,
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
##' @param B a matrix or dataframe of covariates. In an RCT,
##' can include baseline covariates here to control for baseline
##' confounding
##' @param Y outcome variable. Currently should be in [0, 1], but can
##' be continuous.
##' @param a the treatment for which you want to estimate the
##' conditional effect, 1 or 0. 1 corresponds to the ATT, 0
##' corresponds to the ATN or NDE.
##' @param Q.SL.library SuperLearner library for estimating Q
##' @param g.SL.library SuperLearner library for estimating g
##' @param family gaussian for a continuous outcome, binomial for a
##' biniomial outcome
##' @param tol convergence criterion, set to something small
##' @param maxiter maximum number of targetting steps
##' @param target TRUE to preform bias reduction step. FALSE to return
##' simple G computation based estimate
##' @param verbose TRUE for details for each iteration
##' @param Qbound bound Q away from 0 and 1
##' @param gbound bound g away from 0 and 1
##' @param ... aditional parameters to be passed to (both)
##' SuperLearner calls
##' @return adsf
##' @author Sam Lendle \email{lendle@@stat.berkeley.edu}
##' @export
tmle.cte <- function(A, B, Y, a=0, Q.SL.library, g.SL.library, family=gaussian(), tol=1e-10, maxiter=100, target=TRUE, verbose=FALSE, Qbound=c(1e-10, 1-1e-10), gbound=c(1e-10, 1-1e-10), ...) {

  Aa <- as.numeric(A==a)

  Q.init.fit <- SuperLearner(Y, data.frame(A, B), SL.library=Q.SL.library, family=family, ...)
  
  #Q.init.fit <- glm(Qform, gaussian, dat, subset=(Delta==1))
  #Change to SL fit if you want, but remember that gform only has one variable so some
  #algorithms won't work
  Q.A1 <- predict(Q.init.fit, newdata=data.frame(A=1, B))$fit
  Q.A0 <- predict(Q.init.fit, newdata=data.frame(A=0, B))$fit
  
  
  
  g.init.fit <- SuperLearner(A, B, SL.library=g.SL.library, family=binomial,...)
  g.A1 <- .bound(predict(g.init.fit, newdata=B)$fit, gbound)
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
      Q.up <- suppressWarnings(glm(Y~-1+H1, offset=logit(Q.AA), family=binomial))
      #Q.up <- glm(Y~-1+H1, dat, offset=Q.AA, family=gaussian, subset=(Delta==1))
      g.up <- glm(Aa~-1+H2, offset=logit(g.Aa), family=binomial)
      eps <- c(coef(Q.up), coef(g.up))
      if (is.na(eps[2])) {eps[2] <- 0}

      #Fluctuate predicted outcome and treatment, Q and g
      Q.A1 <- logistic(logit(Q.A1)+eps[1]*H1.A1)
      #Q.A1 <- Q.A1+eps[1]*H1.A1
      Q.A0 <- logistic(logit(Q.A0)+eps[1]*H1.A0)
      #Q.A0 <- Q.A0+eps[1]*H1.A0
      g.Aa <- logistic(logit(g.Aa)+eps[2]*H2)
      
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


regress <- function(Y, X, family=binomial(), method="glm", ...) {
  SL.version <- packageVersion("SuperLearner")$major
  if (method=="glm" || !require(SuperLearner)) {
    if (method=="SL") warning("SuperLearner could not be loaded, using main terms glm", call.=FALSE)
    fit <- glm.fit(X, Y, family=family)
  } else if (method=="SL") {
    if (SL.version==1) {
      fit <- SuperLearner(Y, data.frame(X), family=family, ...)
    }
    else {
      stop("regress hasn't been written for new SL")
    }
  }
  res <- list(fit=fit, method=method, SL.version=SL.version)
  class(res) <- "regress"
  return(res)
}

predict.regress <- function(object, newdata=NULL, ...) {
  if (object$method=="glm") {
    pred <- predict.glm(object$fit, newdata=newdata, type=response)
  }
  else if (object$method=="SL") {
    if (object$SL.version==1) {
      if (is.null(newdata)) {
        pred <- predict(object$fit)$fit
      } else {
        pred <-  predict(object$fit, newdata=newdata)
      }
    }
    else {
      warning("This code (predict for new SL) has not been tested")
      pred <- predict(object$fit, newdata=newdata)
    }
  }
  return(pred)
}
