#' classo_fixed
#'
#' cglasso when the covariate is fixed
#' @param X input1
#' @param Y input2
#' @param lambda input3
#' @param weights input4
#' @param b.init input5
#' @return output
#' @examples
#' example <- classo_fixed(X, Y, lambda, weights)
#'
#' @export
classo_fixed <- function(X, Y, lambda, weights, b.init = NULL){

  n <- nrow(X)
  p <- ncol(X)

  ReX <- Re(X); ImX <- Im(X)
  ReY <- Re(Y); ImY <- Im(Y)

  XX <- rbind(cbind(ReX, -ImX),
              cbind(ImX, ReX))
  YY <- c(ReY,
          ImY)

  grp_id <- rep(1:p, 2)

  if (is.null(b.init)){
    b.old <- rep(0, 2*p)
  } else {
    b.old <- c(Re(b.init), Im(b.init))
  }

  r.old <- (YY - XX %*% b.old)
  e.old <- mean(r.old^2)

  iter_max <- 1000
  tol <- 1e-4
  for (iter in 1:iter_max){
    if (is.infinite(e.old) | is.na(e.old)){

      b.old <- rep(NA, 2*p)
      b <- b.old[seq(p)]+1i*b.old[-seq(p)]

      return(b)
    }

    if (e.old > tol){

      r <- r.old
      b <- b.old

      for (j in 1:p){
        grp.idx <- c(j, j+p)
        rj <- r + XX[,c(j,j+p)] %*% b[c(j,j+p)]
        X.rj <- t(XX[,c(j,j+p)]) %*% rj/(n)
        thrs_jk <- weights[j]*lambda

        # Generalized soft threshold operator
        b[c(j,j+p)] <- gen_soft_thrs(X.rj,thrs_jk)

        r <- rj - XX[,c(j,j+p)] %*% b[c(j,j+p)]
      }
      b.old <- b
      r.old <- (YY - XX %*% b.old)
      e.new <- mean(r.old^2)

      if (is.infinite(e.new)|is.infinite(e.old)|is.na(e.new)|is.na(e.old)){

        b.old <- rep(NA, 2*p)
        b <- b.old[seq(p)]+1i*b.old[-seq(p)]

        return(b)
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

  b <- b.old[seq(p)]+1i*b.old[-seq(p)]
  return(b)
}
