#' Simulated data for the cxreg vignette
#'
#' Simple simulated data, used to demonstrate the features of cxreg
#'
#' These datasets are artificial, and are used to test out some of the
#' features of cxreg.
#' @name beta_CVX
#' @aliases x y beta_CVX
#' @format Data objects used to demonstrate features in the cxreg vignette
#' @keywords datasets
#' @useDynLib cxreg
#' @import methods
#' @import Matrix
#' @import foreach
#' @importFrom utils packageDescription
#' @importFrom graphics abline axis matplot points segments text par plot
#' @importFrom stats approx coef median predict runif weighted.mean family rnorm gaussian
#' @importFrom grDevices rainbow
#' @importFrom Rcpp sourceCpp
#' @importFrom fields image.plot
#' @examples
#'
#' data(QuickStartExample)
#' x <- QuickStartExample$x; y <- QuickStartExample$y
#' classo(x, y)
#'
NULL


#' Internal classo functions
#'
#' @description
#' These are not intended for use by users. \code{lambda.interp} does linear
#' interpolation of the lambdas to obtain a prediction at a new point s.
#' \code{nonzeroCoef} determines in an efficient manner which variables are
#' nonzero in each fit.
#' \code{jerr} (not currently available) prints out error messages from the C++ routines.
#' \code{plotCoef} is called by the \code{plot} method for \code{cxreg}
#' objects. \code{check_dots} is used in \code{coef} and \code{predict} with
#' argument \code{exact=TRUE}, to make sure user supplies original data used to
#' fit the \code{"classo"} object.
#'
#' @name classo-internal
#' @aliases cvtype cvstats#'
#' @author Younghoon Kim
#' @keywords internal
NULL

#' Complex-valued Lasso model paths
#'
#' This package fits complex-valued Lassofor regression using coordinate descent. The algorithm is extremely fast, and exploits sparsity in the input x matrix where it exists.
#' A variety of predictions can be made from the fitted models.
#'
#' \tabular{ll}{ Package: \tab cxreg \cr Type: \tab Package\cr Version: \tab
#' 1.0\cr Date: \tab 2025-04-11 \cr License: \tab What license is it under?\cr }
#' Very simple to use. Accepts \code{x,y} data for regression models, and
#' produces the regularization path over a grid of values for the tuning
#' parameter \code{lambda}. Only 5 functions: \code{classo}\cr
#' \code{predict.classo}\cr \code{plot.classo}\cr
#' \code{coef.classo}
#'
#' @name cxreg-package
#' @author Younghoon Kim, Navonil Deb, Sumanta Basu \cr Maintainer:
#' Younghoon Kim <yk748@cornell.edu>
#' @references Deb, N., Kuceyeski, A., Basu, S. (2024)
#' \emph{Regularized Estimation of Sparse Spectral Precision Matrices},
#' \url{https://arxiv.org/abs/2401}
#' @keywords models regression package
#' @examples
#'
#' x <- array(rnorm(100*20), c(100,20)) + (1+1i) * array(rnorm(100*20), c(100,20))
#' for (j in 1:20) x[,j] <- x[,j] / sqrt(mean(Mod(x[,j])^2))
#' e <- rnorm(100) + (1+1i) * rnorm(100)
#' b <- c(1, -1, rep(0, 18)) + (1+1i) * c(-0.5, 2, rep(0, 18))
#' y <- x %*% b + e
#' fit <- glmnet(x, y)
#' predict(fit, newx = x[1:5, ], s = c(0.01, 0.005))
#' predict(fit, type = "coef")
#' plot(fit, xvar = "lambda")
#'
NULL
