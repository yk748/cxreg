#' Cross-validation for classo
#'
#' Does k-fold cross-validation for classo, produces a plot, and returns a
#' value for \code{lambda}
#'
#' The function runs \code{classo} \code{nfolds}+1 times; the first to get the
#' \code{lambda} sequence, and then the remainder to compute the fit with each
#' of the folds omitted. The error is accumulated, and the average error and
#' standard deviation over the folds is computed.
#'
#' Note that the results of \code{cv.classo} are random, since the folds
#' are selected at random. Users can reduce this randomness by running
#' \code{cv.classo} many times, and averaging the error curves.
#'
#' @param x \code{x} matrix as in \code{classo}.
#' @param y response \code{y} as in \code{classo}.
#' @param weights Observation weights; defaults to 1 per observation
#' @param lambda Optional user-supplied lambda sequence; default is \code{NULL},
#' and \code{classo} chooses its own sequence. Note that this is done for the full model (master sequence), and separately for each fold.
#' The fits are then aligned using the master sequence (see the \code{alignment}
#' argument for additional details). Adapting \code{lambda} for each fold
#' leads to better convergence. When \code{lambda} is supplied, the same sequence
#' is used everywhere.
#' @param nfolds number of folds - default is 10. Although \code{nfolds} can be
#' as large as the sample size (leave-one-out CV), it is not recommended for
#' large dataset. Smallest value allowable is \code{nfolds=3}
#' @param foldid an optional vector of values between 1 and \code{nfolds}
#' identifying what fold each observation is in. If supplied, \code{nfolds} can
#' be missing.
#' @param alignment This is an experimental argument, designed to fix the
#' problems users were having with CV, with possible values \code{"lambda"}
#' (the default) else \code{"fraction"}. With \code{"lambda"} the \code{lambda}
#' values from the master fit (on all the data) are used to line up the
#' predictions from each of the folds. In some cases this can give strange
#' values, since the effective \code{lambda} values in each fold could be quite
#' different. With \code{"fraction"} we line up the predictions in each fold
#' according to the fraction of progress along the regularization. If in the
#' call a \code{lambda} argument is also provided, \code{alignment="fraction"}
#' is ignored (with a warning).
#' @param keep If \code{keep=TRUE}, a \emph{prevalidated} array is returned
#' containing fitted values for each observation and each value of \code{lambda}.
#' This means these fits are computed with this observation and the rest of its fold omitted.
#' The \code{foldid} vector is also returned. Default is keep=FALSE.
#' @param parallel If \code{TRUE}, use parallel \code{foreach} to fit each
#' fold. Must register parallel before hand, such as \code{doMC} or others.
#' Currently it is unavailable.
#' @param trace.it If \code{trace.it=1}, then progress bars are displayed;
#' useful for big models that take a long time to fit. Limited tracing if
#' \code{parallel=TRUE}
#' @param \dots Other arguments that can be passed to \code{classo}
#'
#' @return an object of class \code{"cv.classo"} is returned, which is a list
#' with the ingredients of the cross-validation fit.
#' \item{lambda}{the values of \code{lambda} used in the fits.}
#' \item{cvm}{The mean cross-validated error - a vector of length \code{length(lambda)}.}
#' \item{cvsd}{estimate of standard error of \code{cvm}.}
#' \item{cvup}{upper curve = \code{cvm+cvsd}.}
#' \item{cvlo}{lower curve = \code{cvm-cvsd}.}
#' \item{nzero}{number of non-zero coefficients at each \code{lambda}.}
#' \item{name}{a text string indicating type of measure for plotting purposes).}
#' \item{classo.fit}{a fitted classo object for the full data.}
#' \item{lambda.min}{value of \code{lambda} that gives minimum \code{cvm}.}
#' \item{lambda.1se}{largest value of \code{lambda} such that error is within 1 standard error of the minimum.}
#' \item{fit.preval}{if \code{keep=TRUE}, this is the array of pre-validated fits. Some entries can be \code{NA},
#' if that and subsequent values of \code{lambda} are not reached for that fold}
#' \item{foldid}{if \code{keep=TRUE}, the fold assignments used}
#' \item{index}{a one column matrix with the indices of \code{lambda.min} and \code{lambda.1se} in the sequence of coefficients, fits etc.}
#' @author Navonil Deb, Younghoon Kim, Sumanta Basu \cr Maintainer: Younghoon Kim
#' \email{yk748@cornell.edu}
#' @seealso \code{classo} and \code{plot} and \code{coef} methods for \code{"cv.classo"}.
#' @examples
#'
# set.seed(1010)
# n = 1000
# p = 200
# x = array(rnorm(n*p), c(n,p)) + (1+1i) * array(rnorm(n*p), c(n,p))
# for (j in 1:p) x[,j] = x[,j] / sqrt(mean(Mod(x[,j])^2))
# e = rnorm(n) + (1+1i) * rnorm(n)
# b = c(1, -1, rep(0, p-2)) + (1+1i) * c(-0.5, 2, rep(0, p-2))
# y = x %*% b + e
#' cv.test = cv.classo(x,y)
#'
#' @export cv.classo
cv.classo <- function (x, y,
                       weights=NULL,
                       lambda = NULL,
                       nfolds = 10,
                       foldid=NULL,
                       type.measure="mse",
                       alignment=c("lambda","fraction"),
                       keep = FALSE,
                       parallel = FALSE,
                       trace.it=0, ...){

  # ------------------------------------------------ #
  type.measure <- match.arg(type.measure)
  alignment <- match.arg(alignment)
  if (!is.null(lambda) && length(lambda) < 2){
    stop("Need more than one value of lambda for cv.glasso")
  }
  if (!is.null(lambda) && alignment=="fraction"){
    warning("fraction of path alignment not available if lambda given as argument; switched to alignment=`lambda`")
    alignment <- "lambda"
  }

  # ------------------------------------------------ #
  N <- nrow(x)
  if (is.null(weights)){
    weights <- rep(1, N)
  } else {
    weights <- as.double(weights)
  }

  # ------------------------------------------------ #
  y <- drop(y)
  cv.call <- classo.call <- match.call(expand.dots = TRUE)
  which <- match(c("type.measure","nfolds","foldid","keep"), names(classo.call), FALSE)
  if (any(which)){
    classo.call <- classo.call[-which]
  }

  # ------------------------------------------------ #
  classo.call[[1]] <- as.name("classo")
  if(classo.control()$itrace){
    trace.it <- 1
  }else{
    if(isTRUE(trace.it == 1)){
      classo.control(itrace=1)
      on.exit(classo.control(itrace=0))
    }
  }

  # ------------------------------------------------ #
  if (is.null(foldid)){
    foldid <- sample(rep(seq(nfolds), length = N))
  }else {
    nfolds <- max(foldid)
  }

  if (nfolds < 3){
    stop("nfolds must be bigger than 2; nfolds=10 recommended")
  }

  # ------------------------------------------------ #
  # Call cv.classo.raw
  cv.classo.raw(x,y,weights,lambda,type.measure,nfolds,foldid,
                alignment,keep,parallel,trace.it,classo.call,cv.call)
}
