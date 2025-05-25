#' fit a complex-valued graphical lasso
#'
#' Fit a complex-valued graphical lasso for spectral precision matrix (inverse spectral density matrix) via a complex variable-wise coordinate descent algorithm for classo with covariates update.
#'
#' The sequence of models implied by \code{lambda} is fit by coordinate descent.
#'
#' @param S p x p-dimensional symmetric spectral density (or spectral coherence) matrix. S is considered as being computed by average smoothed periodogram (the bandwidth is computed by using the given nobs).
#' @param scale The provided S is spectral density (covariance) or spectral coherence (coherence). Default is \code{covariance}.
#' @param D The p x p-dimensional diagonal matrix with spectral densities as the diagonal entries. Default is \code{NULL}. Do not provide this input if scale is \code{covariance} and type is \code{I}; the algorithm will stop. If not provided when scale is \code{TRUE}, identity matrix is chosen.
#' @param type Logical flag to choose the formulation to solve. Default is \code{I}. If type is \code{I}, the algorithm solves CGLASSO-I in the reference, 
#' \deqn{D^{-1/2}(\arg\min_{\Theta}Tr[\hat{R}\hat{\Theta}]-\log\det\Theta + \sum_{i\neq j}|\Theta_{ij}|)D^{-1/2}} for the given $D$. If type is \code{II}, the algorithm solves CGLASSO-II in the reference. It is for each iterative classo with covariate update, the scale matrix D is multiplied. Please see the reference for the details of the iterative updates.
#' @param nobs Number of observations used in computation of the spectral density matrix S. This quantity is need to compute the Fourier frequency, extended BIC, and bandwidth for the average smoothed periodogram. 
#' @param lambda A user supplied \code{lambda} sequence.
#' Typical usage is to have the program compute its own \code{lambda} sequence based on
#' \code{nlambda} and \code{lambda.min.ratio}.
#' Supplying a value of \code{lambda} overrides this.
#' WARNING: use with care. Avoid supplying a single value for \code{lambda}
#' @param nlambda The number of \code{lambda} values - default is 50.
#' @param lambda.min.ratio Smallest value for \code{lambda}, as a fraction of
#' \code{lambda.max}, the (data derived) entry value (i.e. the smallest value
#' for which all coefficients are zero). The default depends on the sample size
#' \code{nobs} relative to the number of variables \code{nvars}.
#' If \code{nobs > p}, the default is \code{0.0001}, close to zero.
#' If \code{nobs < p}, the default is \code{0.01}.
#' @param W.init Logical flag whether the initially estimated spectral density matrix is given. Default is \code{NULL}.
#' @param stopping_rule Logical flag if the algorithm is terminated by stopping rule. If the algorithm is early terminated ,not all estimates for initially designated lambdas are explored. 
#' @param stop_criterion Stopping criterion for early termination. Default is \code{EBIC} (Extended BIC). Alternatively, \code{AIC} (AIC) and \code{RMSE} (root mean squared error between two consecutive estimates) can be used.
#' @param maxit Maximum number of iterations of both outer and inner loops. Default 500.
#' @param thresh Convergence threshold for coordinate descent. Default is 1e-4.
#' @param trace.it If \code{trace.it=1}, then a progress bar is displayed;
#' useful for big models that take a long time to fit.
#' 
#' @return An object with class "cglassofit" and "cglasso".
#' \item{stop_arr}{Sequence of values of information criterion for a fixed lambda.}
#' \item{stop_criterion}{Stopping criterion used.}
#' \item{min_index}{The index for lambda that minimizes the value of the information criterion.}
#' \item{lambda_grid}{Sequence of lambdas used.}
#' \item{Theta_hat}{Estimated inverse spectral matrix for each fixed lambda. It is provided in the list.}
#' \item{type}{Type of the formulation used, either CGALSSO-I or CGLASSO-II.}
#' \item{scale}{Whether the spectral density matrix (covariance) or spectral coherence (coherence) is given.}
#' \item{D}{Used scale diagonal matrix.}
#' 
#' @author Navonil Deb, Younghoon Kim, Sumanta Basu \cr Maintainer: Younghoon Kim
#' \email{yk748@cornell.edu}
#' @references Deb, N., Kuceyeski, A., Basu, S. (2024)
#' \emph{Regularized Estimation of Sparse Spectral Precision Matrices},
#' \url{https://arxiv.org/abs/2401.11128}.
#' @keywords models complex-valued precision matrix
#' @examples
#' p <- 30
#' n <- 500
#' C <- diag(0.7, n)
#' C[row(C) == col(C) + 1] <- 0.3  
#' C[row(C) == col(C) - 1] <- 0.3  
#' Sigma <- solve(C)
#' set.seed(1010)
#' X_t <- rmvnorm(n = n, mean = rep(0, p), sigma = Sigma)
#' d_j <- dft.X(X_t,j,m)
#' f_j_hat <- t(d_j) %*% Conj(d_j) / (2*m+1)
#' fit <- cglasso(S=f_j_hat, scale="covariance", nobs=n)
#' 
#' #' @export cglasso
cglasso <- function(S,
                    scale = "covariance",
                    D = NULL,
                    type = "I",
                    nobs,
                    lambda = NULL,
                    nlambda = 50,
                    lambda.min.ratio = ifelse(nobs<p,1e-2,1e-4),
                    W.init = NULL,
                    stopping_rule = TRUE,
                    stop_criterion = c("EBIC","AIC","RMSE"),
                    maxit = 500,
                    thresh = 1e-4,
                    trace.it = 0,...){
  
  this.call <- match.call()
  ####################################################################
  # check for NAs in x
  if(any(is.na(S))){
    stop("S has missing values; S needs to be complete")
  }else{
    p1 <- dim(S)[1]; p2 <- dim(S)[2]
    if (p1 != p2){
      stop("S must be square")
    }else{
      p <- p1
    }
  }
  # check for NAs in nobs
  if(is.null(nobs)){
    stop("nobs is missing; the number of observations needs to be provided")
  }
  # check if the D is not given when scale is covariance & type is I 
  # or if D is properly given under the different settings.
  if (scale=="covariance" & type=="I" & !is.null(D)){
    stop("D should not be provided if type I formulation is considered with covariance scale.")
  }else if (!any(is.null(D))){
    if ( (dim(D)[1] != p) | (dim(D)[2] != p) ){
      stop("D must be p by p matrix")
    }else if (!all(D[!diag(nrow(p))] == 0)){
      stop("D must be diagonal")
    }
  }else{
    D <- diag(1,p)
  }
  # Check stopping rule & stopping criterion
  if (is.null(stop_criterion)){
    stop_criterion <- "EBIC"
  }
  
  
  ####################################################################
  fit_cglasso <- cglasso.path(S,
                              scale,
                              D,
                              type,
                              nobs,
                             lambda,
                             nlambda,
                             lambda.min.ratio,
                             W.init,
                             scale,
                             stopping_rule,
                             stop_criterion,
                             maxit = maxit,
                             thresh = thresh,
                             trace.it = trace.it)
  
  fit_cglasso$call <- this.call
  class(fit_cglasso) <- "cglasso"
  fit_cglasso
}
