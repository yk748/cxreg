#' fit a complex-valued lasso
#'
#' Fit a complex-valued lasso formulation via complex update coordinate descent algorithm.
#' By defining a field isomophism between complex values and its 2 by 2 representation,
#' it enables to update each coordinate of the solution as a regular coordinate descent algorithm.
#'
#' The sequence of models implied by \code{lambda} is fit by coordinate descent.
#' For \code{family="gaussian"} this is the lasso sequence.
#'
#' @param x input matrix, of dimension nobs x nvars; each row is an observation vector.
#' Requirement: \code{nvars >1}; in other words, \code{x} should have 2 or more columns.
#' @param y response variable.
#' @param family To match the conventional functions used in \code{glmnet}, "gaussian" is assumed. This has to be eliminated.
#' @param weights observation weights. Default is 1 for each observation.
#' @param nlambda The number of \code{lambda} values - default is 100.
#' @param lambda.min.ratio Smallest value for \code{lambda}, as a fraction of
#' \code{lambda.max}, the (data derived) entry value (i.e. the smallest value
#' for which all coefficients are zero). The default depends on the sample size
#' \code{nobs} relative to the number of variables \code{nvars}.
#' If \code{nobs > nvars}, the default is \code{0.0001}, close to zero.
#' If \code{nobs < nvars}, the default is \code{0.01}.
#' A very small value of \code{lambda.min.ratio} will lead to a saturated fit
#' in the \code{nobs < nvars} case.
#' @param lambda A user supplied \code{lambda} sequence.
#' Typical usage is to have the program compute its own \code{lambda} sequence based on
#' \code{nlambda} and \code{lambda.min.ratio}.
#' Supplying a value of \code{lambda} overrides this.
#' WARNING: use with care. Avoid supplying a single value for \code{lambda}
#' (for predictions after CV use \code{predict()} instead).
#' Supply instead a decreasing sequence of \code{lambda} values.
#' \code{classo} relies on its warms starts for speed, and its often faster to fit a whole path than compute a single fit.
#' @param standardize Logical flag for x variable standardization, prior to
#' fitting the model sequence. The coefficients are always returned on the
#' original scale. Default is \code{standardize=FALSE}.
#' @param intercept Should intercept(s) set to zero (default=FALSE) or be fitted (TRUE).
#' @param thresh Convergence threshold for coordinate descent.
#' Each inner coordinate-descent loop continues until the maximum change in the objective
#' after any coefficient update is less than \code{thresh} times the null
#' deviance. Defaults value is \code{1E-7}.
#' @param trace.it If \code{trace.it=1}, then a progress bar is displayed;
#' useful for big models that take a long time to fit.
#'
#' @author Navonil Deb, Younghoon Kim, Sumanta Basu \cr Maintainer: Younghoon Kim
#' \email{yk748@cornell.edu}
#' @references Deb, N., Kuceyeski, A., Basu, S. (2024)
#' \emph{Regularized Estimation of Sparse Spectral Precision Matrices},
#' \url{https://arxiv.org/abs/2401.11128}.
#' @keywords models complex-valued regression
#' @examples
#' set.seed(1010)
#' n = 1000
#' p = 200
#' x = array(rnorm(n*p), c(n,p)) + (1+1i) * array(rnorm(n*p), c(n,p))
#' for (j in 1:p) x[,j] = x[,j] / sqrt(mean(Mod(x[,j])^2))
#' e = rnorm(n) + (1+1i) * rnorm(n)
#' b = c(1, -1, rep(0, p-2)) + (1+1i) * c(-0.5, 2, rep(0, p-2))
#' y = x %*% b + e
#' fit.test = classo(x, y)
#'
#' @export classo
classo <- function(x,y,
                   weights=NULL,
                   family=gaussian(),
                   lambda=NULL,
                   nlambda=100,
                   lambda.min.ratio=ifelse(nobs<nvars,1e-2,1e-4),
                   standardized=FALSE,
                   intercept=FALSE,
                   maxit=100000,
                   trace.it=0,...){

  # ------------------------------------------------ #
  # x, # included
  # y, # included
  # family="gaussian", # included
  # weights=NULL, # included
  # offset=NULL, # deleted
  # alpha=1.0, # deleted
  # nlambda=100, # included
  # lambda.min.ratio=ifelse(nobs<nvars,1e-2,1e-4), # included
  # lambda=NULL, # included
  # standardize=TRUE, # included
  # intercept=FALSE, # included
  # thresh=1e-7, # included
  # dfmax=nvars+1, # deleted
  # pmax=min(dfmax*2+20,nvars), # deleted
  # exclude=NULL, # deleted
  # penalty.factor=rep(1,nvars), # deleted
  # lower.limits=-Inf, # deleted
  # upper.limits=Inf, # deleted
  # maxit=100000, # included
  # type.gaussian=ifelse(nvars<500,"covariance","naive"), # deleted
  # type.logistic=c("Newton","modified.Newton"), # deleted
  # standardize.response=FALSE, # deleted
  # type.multinomial=c("ungrouped","grouped"), # deleted
  # relax=FALSE, # deleted
  # trace.it=0 # included
  # ------------------------------------------------ #

  # YK: gaussian() and related functions should be replaced.

  this.call <- match.call()
  # ------------------------------------------------ #
  # Need to do this first so defaults in call can be satisfied
  np <- dim(x)

  # check dims
  if(is.null(np)|(np[2]<=1)){
    stop("x should be a matrix with 2 or more columns")
  }
  nobs <- np[1]
  nvars <- np[2]

  # check for NAs in x
  if(any(is.na(x))){
    stop("x has missing values; consider using makeX() to impute them")
  }

  # check for NAs in weights
  if(is.null(weights)){
    weights <- rep(1,nobs)
  }else if(length(weights)!=nobs){
    stop(paste("number of elements in weights (",length(weights),")
               not equal to the number of rows of x (",nobs,")",sep=""))
  }

  # ------------------------------------------------ #
  # See whether its a call to glmnet or to glmnet.path, based on family arg
  # note: we exclude the rest of condition checking steps.
  # note: here we have only one type of family.
  fit <- classo.path(x,y,family = gaussian(),
                     weights=NULL,
                     standardized=FALSE,
                     lambda=NULL,
                     nlambda=100,
                     lambda.min.ratio=ifelse(nobs<nvars,1e-2,1e-4),
                     intercept = FALSE,
                     thresh = 1e-10,
                     maxit = 100000,
                     trace.it = 0)

  fit$call <- this.call
  fit
}

