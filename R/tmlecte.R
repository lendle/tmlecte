##' Prints results from tmle.cte.
##'
##' 
##' @title \code{print.cte} function
##' @param x an ojbect of class \code{cte}
##' @param ... ignored
##' @return Nothing
##' @author Sam Lendle
##' @method print cte
##' @export
print.cte <- function(x, ...) {
  cat(x$estimand, ": ", x$psi,
      "\nEstimated Variance: ", x$var.psi,
      "\n95% Confidence interval: (", x$CI[1], ", ", x$CI[2],")",
      "\nP-value for H0: psi=0 vs. H1: psi!=0:", format.pval(x$pvalue, eps=0.001), "\n", sep="")
}



##' Estimates the statistical parameter E(E(Y|A=1,W)-E(Y|A=0,W)|A=a) using Targetted Maximum Likelihood.
##'
##' If \code{target} is \code{FALSE}, the targetting step is not done, and a G-computation type estimate based only on the inital estimate of Q is returned.  In this case, the influence curve, variance estimate, confidence interval and p-value are not calculated.
##' 
##' Q.SL.library defaults to ('SL.glm', 'SL.step', 'SL.glm.interaction')
##' g.SL.library Defaults to ('SL.glm', 'SL.step', 'SL.glm.interaction')
##' This choice is simply because these algorithms are included in the base R installation. See SuperLearner help files for further information.
##' @title tmle.cte
##' @param Y continuous or binary outcome variable
##' @param A binary treatment indicator, \code{1} - treatment, \code{0} - control
##' @param B vector, matrix, or dataframe containing covariates to control for. Baseline covariates for the ATT,
##' baseline covariates and mediator for NDEU
##' @param a the treatment for which you want to estimate the
##' conditional effect, \code{1} or \code{0}. \code{1} corresponds to the ATT, \code{0}
##' corresponds to the ATN or NDE.
##' @param Delta indicator of missing outcome or treatment assignment. \code{1} - observed, \code{0} - missing
##' @param Q optional \emph{nx2} matrix of initial values for Q portion of the likelihood, \emph{(E(Y|A=0,B), E(Y|A=1,B))}
##' @param Qform optional regression formula for estimation of \emph{E(Y|A,B)}, suitable for call to \code{glm}
##' @param Q.SL.library optional vector of prediction algorithms to use for \code{SuperLearner} estimation of initial Q
##' @param Qbounds vector of upper and lower bounds on Y and predicted values for initial \code{Q}. Defaults to the range of \code{Y}, widened by 10\% of the min and max values.
##' @param g1W optional vector of conditional treatment assingment probabilities, \emph{P(A=1|B)}
##' @param gform optional regression formula of the form A~W, if specified this overrides the call to \code{SuperLearner}
##' @param gbound value between (0,1) for truncation of predicted probabilities. See \code{Details} section for more information
##' @param pDelta1 optional matrix of conditional probabilities for missingness mechanism with dimension \emph{nx2} \emph{P(Delta=1|A=0,B)}, \emph{P(Delta=1|A=1,B)}
##' @param g.Deltaform optional regression formula of the form \code{Delta~A+B}, if specified this overrides the call to \code{SuperLearner}
##' @param g.SL.library optional vector of prediction algorithms to use for \code{SuperLearner} estimation of \code{g1W} or \code{pDelta1}
##' @param family family specification for working regression models, generally 'gaussian' for continuous outcomes (default), 'binomial' for binary outcomes
##' @param fluctuation 'logistic' (default), or 'linear'
##' @param alpha used to keep predicted initial values bounded away from (0,1) for logistic fluctuation
##' @param id optional subject identifier
##' @param target \code{TRUE} (default) to preform bias reduction (targeting)
##' step. \code{FALSE} to return simple G computation based substitution
##' estimate
##' @param tol convergence criterion, defaults to square root of machine epsilon
##' @param maxiter maximum number of iterations in the targeting step. Default is 100.
##' @param verbose status messages printed if set to \code{TRUE} (default=\code{FALSE})
##' @return An object of class \code{tmle.cte} which is a list with the following components:
##' \describe{
##' \item{\code{psi}}{Parameter estimate}
##' \item{\code{var.psi}}{Influence curve based estimate of the variance of \code{psi}}
##' \item{\code{CI}}{95\% confidence interval for \code{psi}}
##' \item{\code{pvalue}}{p-value from a two-sided test of \code{H0: psi_0=0}}
##' \item{\code{IC}}{A vector with the value of the estimated influence curve for each observation}
##' \item{\code{estimand}}{Character string describing the estimand, depends on the value of \code{a}}
##' \item{\code{Qinit}}{Initial estimate of Q from the \code{estimateQ} function}
##' \item{\code{Qstar}}{Matrix of final estimates of Q after targeting}
##' \item{\code{ginit}}{Initital estimate of g from the \code{estimateg} function}
##' \item{\code{g1Wstar}}{Vector of final estiates of g after targeting}
##' \item{\code{g.Delta}}{Estimate of P(Delta=1|A,W) from the \code{estimateg} function}
##' }
##' If convergence fails, all items are either \code{NA} or \code{NULL}.
##' If \code{target==FALSE}, some items are \code{NA} or \code{NULL}.
##' @author Sam Lendle \email{lendle@@stat.berkeley.edu}
##' @seealso \code{\link{tmle.att}}, \code{\link{tmle.ndeu}}
##' @export
tmle.cte <- function(Y, A, B,
                     a=0,
                     Delta=NULL,
                     Q=NULL,
                     Qform=NULL,
                     Q.SL.library=c("SL.glm", "SL.step", "SL.glm.interaction"),
                     Qbounds=NULL,
                     g1W=NULL,
                     gform=NULL, 
                     gbound=0.001,
                     pDelta1=NULL,
                     g.Deltaform,
                     g.SL.library=c("SL.glm", "SL.step", "SL.glm.interaction"),
                     family=gaussian(),
                     fluctuation=c("logistic", "linear"),                     
                     alpha=0.995,
                     id=1:length(Y),
                     target=TRUE,
                     tol=sqrt(.Machine$double.eps),
                     maxiter=100,
                     verbose=FALSE
                     ) {

  ## #####################
  ## Checking arguments ##
  ## #####################
  
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

  if ((is.null(Delta) || all(Delta==1)) && all(!is.na(Y))) {
    missing.outcome <- FALSE
    Delta <- rep(1, length(Y))
  } else {
    missing.outcome <- TRUE
    if (is.null(Delta)) {
      warning("NA values found in Y, but Delta is NULL.  Setting Delta=0 if Y is NA, and 1 otherwise")
      Delta <- as.numeric(!is.na(Y))
    }
    if (length(Delta) != length(Y)) {
      warning("Delta is the wrong length. Setting Delta=0 if Y is NA, and 0 otherwise")
      Delta <- as.numeric(!is.na(Y))
    }
    if (any(is.na(Y)&(Delta==1))) {
      warning("Some values of Y are NA where Delta is 1.  Setting those values of Delta to 0")
      Delta[is.na(Y)] <- 0
    }
    Y[Delta==0] <- NA    
  }

  if (!is.null(Q) && (nrow(Q) != length(Y) || ncol(Q) != 2)) {
    stop ("User specified Q should have dimension length(Y) x 2")
  }
  if (!is.null(g1W) && length(g1W) != length(Y)) {
    stop ("The length of user specified g1W should be the same as the length of Y")
  }
  if (missing.outcome && !is.null(pDelta1) && (nrow(pDelta1) != length(Y) || ncol(pDelta1) != 2)) {
    stop ("User specified pDelta1 should have dimension length(Y) x 2")
  }

  fluctuation <- match.arg(fluctuation)

  ## ##########################
  ## Done checking arguments ##
  ## ##########################  

  ## ###############################
  ## Set up and initial estimates ##
  ## ###############################
  
  colnames(B) <- .setColnames(colnames(B), NCOL(B), "B")
  
  ## Aa is a dummy treatment variable, A if a=1, 1-A if a=0
  ## for use in the update of g, when fitting P(A=a | W)
  Aa <- as.numeric(A==a)

  maptoYstar <- fluctuation=="logistic"

  if (length(gbound)==1) gbound <- c(gbound, 1-gbound)   

  stage1 <- .initStage1(Y, A, Q, Q.Z1=NULL, Delta, Qbounds, alpha, maptoYstar, family$family)		
  
  Q <- estimateQ(Y=stage1$Ystar,
                 Z=0, 
                 A=A,
                 W=B,
                 Delta=Delta,
                 Q=stage1$Q,
                 Qbounds=stage1$Qbounds,
                 Qform=Qform,
                 maptoYstar=maptoYstar, 
                 SL.library=Q.SL.library,
                 cvQinit=FALSE,
                 family=family$family,
                 id=id,
                 verbose=verbose)

  g <- estimateG(d=data.frame(A, B),
                 g1W=g1W,
                 gform=gform,
                 SL.library=g.SL.library,
                 id=id,
                 verbose=verbose,
                 message="treatment mechanism",
                 outcome="A")

  g.Delta <- estimateG(d=data.frame(Delta, Z=0, A, B),
                       g1W=pDelta1,
                       gform=g.Deltaform,
                       SL.library=g.SL.library,
                       id=id,
                       verbose=verbose,
                       message="missingness mechanism",
                       outcome="D")


  ## Flag denoteing if convergence fails
  fail <- FALSE
  Qstar <- Q$Q
  
  if (target) {
    iter <- 0
    done <- FALSE
    prev.crit <- -1
    ## gaW is P(A=a | W)
    gaW <- a*g$g1W + (1-a)*(1-g$g1W)
    if (Q$family=="binomial") {
      Qupfam <- "quasibinomial"
    } else {
      Qupfam <- Q$family
    }
  } else {
    ## to skip targeting step
    done <- TRUE
  }

  ## ###################################
  ## Done setup and initial estimates ##
  ## ###################################

  ## ##################
  ## Targeting step ##
  ## ##################

  while(!done) {
    iter=iter+1
    
    ##Calculate target parameter at current step
    if (maptoYstar) {
      psi <- mean(plogis(Qstar[A==a, "Q1W"]) - plogis(Qstar[A==a, "Q0W"]))
    } else {
      psi <- mean(Qstar[A==a, "Q1W"] - Qstar[A==a, "Q0W"])      
    }
    
    ##bound g away from 0 and 1 so logit(g) can be calculated, and clever covars
    gaW <- .bound(gaW, gbound)

    if (a==1) {
      g1W <- gaW
      g0W <- 1-g1W
    } else {
      g0W <- gaW
      g1W <- 1-g0W
    }
    
    ##Calculate covariates for fluctuation
    H1.1W <- (1/g.Delta$g1W[,"Z0A1"])*(gaW/g1W)
    H1.0W <- (1/g.Delta$g1W[,"Z0A0"])*(-gaW/g0W)
    H1 <- A*H1.1W + (1-A)*H1.0W
    H2 <- (Qstar[,"Q1W"]-Qstar[,"Q0W"] - psi)


    ##Estimate epsilon for logistic fluctuation of Q and g
    Q.up <- (glm(stage1$Ystar~-1+H1, offset=Qstar[,"QAW"], family=Qupfam))
    g.up <- glm(Aa~-1+H2, offset=qlogis(gaW), family=binomial)
    eps <- c(coef(Q.up), coef(g.up))
    if (is.na(eps[2])) {eps[2] <- 0}

    ##Fluctuate predicted outcome and treatment, Q and g
    Qstar[,"Q1W"] <- Qstar[,"Q1W"]+eps[1]*H1.1W
    Qstar[,"Q0W"] <- Qstar[,"Q0W"]+eps[1]*H1.0W
    Qstar[,"QAW"] <- A*Qstar[,"Q1W"] + (1-A)*Qstar[,"Q0W"]
    if (fluctuation=="logistic") {
      Qstar <- .bound(Qstar, qlogis(stage1$Qbounds))
    } else {
      Qstar <- .bound(Qstar, stage1$Qbounds)
    }
   
    gaW <- plogis(qlogis(gaW)+eps[2]*H2)
    
    crit <- max(abs(na.omit(eps)))
    if (crit <= tol || iter>=maxiter) {
      done=TRUE
    }

    ## If not done and critical value has not improved (possibly because all values in Qstar are
    ## at the boundary) then convergence failed
    ## Note that if Q is not bounded, the offset can be +/-Inf so convergence fails either way...
    if(done==FALSE && crit==prev.crit) {
      done=TRUE
      fail=TRUE
      warning("Convergence criterion failed to improve between iterations")
    }
    prev.crit <- crit
    
    if (verbose) cat("iter:", iter, "eps:", eps, "crit:", crit, "psi:", psi, "\n")
  }

  if (maptoYstar) {
    Qstar <- plogis(Qstar)*diff(stage1$ab)+stage1$ab[1] 
    Q$Q <- plogis(Q$Q)*diff(stage1$ab)+stage1$ab[1]
  }

  if (fail) {
    res <- list(psi=NA,
                var.psi=NA,
                CI=NA,
                pvalue=NA,
                IC=NA,
                estimand=NULL,
                Qinit=NULL,
                Qstar=NULL,
                ginit=NULL,
                g1Wstar=NULL,
                g.Delta=NULL
                )
  }

  if (!fail) {
    #Calculate final parameter of interest after convergence
    psi <-  mean(Qstar[A==a, "Q1W"] - Qstar[A==a, "Q0W"])
    if (target) {
      p.Aa <- mean(Aa)
    #To simplify calculation of IC
      Y[Delta==0] <- 0
      IC <- ((Delta)*(A*gaW/(g.Delta$g1W[,"Z0A1"]*g1W) - (1-A)*gaW/(g.Delta$g1W[,"Z0A0"]*g0W))*(Y-Qstar[, "QAW"])+Aa*(Qstar[,"Q1W"]-Qstar[,"Q0W"] - psi))/p.Aa
      if (verbose) {
        cat("Mean of the influence curve:", mean(IC), "\n")
      }
      var.psi <- var(IC)/length(Y)
      CI <- psi + c(-1.96, 1.96)*sqrt(var.psi)
      pvalue <- 2*(1-pnorm(abs(psi/sqrt(var.psi))))
      if (a==1) {
        g1Wstar <- gaW
      } else {
        g1Wstar <- 1-gaW
      }
    } else {
      ## If targeting is skipped, do not calculate variance, CI, p-value or IC
      var.psi <- IC <- CI <- pvalue <- NA
      Qstar <- g1Wstar <- NULL
    }
    estimand <- paste("E(E(Y|A=1,B)-E(Y|A=0,B)|A=", a, ")", sep="")
    res <- list(psi=psi,
                var.psi=var.psi,
                CI=CI,
                pvalue=pvalue,
                IC=IC,
                estimand=estimand,
                Qinit=Q,
                Qstar=Qstar,
                ginit=g,
                g1Wstar=g1Wstar,
                g.Delta)
  }
  class(res) <- c("cte")
  res
}
