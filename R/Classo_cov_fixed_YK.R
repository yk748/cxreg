#' classo_cov_fixed
#'
#' cglasso with covariate update
#' @param X input1
#' @param Y input2
#' @param lambda input3
#' @param weights input4
#' @param b.init input5
#' @return output
#' @examples
#' example <- classo_cov_fixed(X, Y, lambda, weights)
#'
#' @export
classo_cov_fixed <- function(X, Y, lambda, weights, b.init = NULL){

  n <- nrow(X)
  p <- ncol(X)

  XX <- t(Conj(X)) %*% X /n
  XY <- t(Conj(X)) %*% Y /n

  grp_id <- rep(1:p, 2)

  if (is.null(b.init)){
    b.old <- complex(p, real = 0, imaginary = 0)
  } else {
    b.old <- b.init
  }

  X.r.old <- XY - XX %*% b.old
  e.old <- mean(Mod(X.r.old)^2)

  iter_max <- 1000
  tol <- 1e-4
  for (iter in 1:iter_max){

    if (is.infinite(e.old) | is.na(e.old)){

      b.old <- complex(p, real = NA, imaginary = NA)
      return(b.old)
    }

    if (e.old > tol){
      X.r <- X.r.old
      b <- b.old

      for (j in 1:p){

        X.rj <- X.r + XX[,j]*b[j]
        Xj.rj <- X.rj[j]
        thrs_jk <- weights[j]*lambda

        # Soft threshold operator (typical lasso operator for a single entry)
        b[j] <- soft_thrs_single(Xj.rj,XX[,j],thrs_jk,scale)

        X.r <- X.rj - XX[,j] * b[j]

      }
      b.old <- b
      X.r.old <- X.r
      e.new <- mean(Mod(X.r.old)^2)

      if (is.infinite(e.new)|is.infinite(e.old)|is.na(e.new)|is.na(e.old)){

        b.old <- complex(p, real = NA, imaginary = NA)
        return(b.old)
      }

      if (abs(e.new - e.old) < tol){
        break
      } else {
        e.old <- e.new
      }

    } else {
      break
    }
  }

  return(b.old)
}
