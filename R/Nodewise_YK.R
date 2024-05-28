#' classo_warm
#'
#' Assume columns of X have norm = sqrt(N)
#' @param Xt input1
#' @param j input2
#' @param method input3
#' @param type input4
#' @param scale input5
#' @param weights input6
#' @param pathwise input7
#' @return output
#' @examples
#' example <- nodewise_YK(Xt, j, method = "or", type = "ols", scale = TRUE, weights=NULL, pathwise=FALSE)
#'
#' @export
nodewise_YK <- function(Xt, j, method = "or", type = "ols", scale = TRUE, weights=NULL, pathwise=FALSE){

  n <- dim(Xt)[1]
  p <- dim(Xt)[2]

  # Construct \mathscf{Z}:
  m <- round(sqrt(n))
  N <- 2*m+1
  Z_scf <- dft.X(Xt, j, m) # 2m+1 x p

  # Columnwise scaling:
  if (scale == TRUE){
    Z_scale <- Z_scf
    Z_scale <- mapply(x=1:p,
                      function(x) Z_scf[,x] <- (Z_scf[,x]-mean(Z_scf[,x])))
    Z_scale <- mapply(x=1:p,
                      function(x) Z_scale[,x] <- sqrt(N)*Z_scale[,x]/norm(Z_scale[,x],"2"))

    G_container <- array(NA, dim=c(p,p))
    diag(G_container) <- 1

  }else{
    # If scale == FALSE,
    # Weight must not be NULL,
    # And Class_fixed & Class_cov_fixed must take scale as an input,
    # Thus Complex operators must be modified as well.
    Z_scale <- Z_scf
    f_hat_avg <- t(Z_scale) %*% Z_scale / (2*m+1)

    det <- prod(eigen(f_hat_avg)$values)
    G_container <- array(NA, dim=c(p,p))
    diag(G_container) <- diag(f_hat_avg)/det
  }

  # Weights:
  if (is.null(weights)){
    W <- matrix(1,nrow=p,ncol=p)
  }else{
    W <- weights
  }

  # Ready:
  eps <- 1e-3; Lambda_length <- 50
  b_hat <- array(NA,dim=c(p,(p-1)))

  mse <- array(NA, dim=c(Lambda_length,p))
  bic_beta <- array(NA, dim=c(Lambda_length,p))
  for(k in 1:p){

    yk <- Z_scale[ ,k]
    Xk <- Z_scale[ ,-k]

    Lambda_max <- max(abs(t(Xk) %*% yk))/N
    Lambda_min <- eps*Lambda_max
    Lambdas <- exp(seq(log(Lambda_min),log(Lambda_max),length=Lambda_length))

    vec <- vector("numeric",length=p)
    arr_beta <- matrix(NA, nrow=length(Lambdas), ncol=(p-1))
    for(l in 1:length(Lambdas)){

      if (pathwise == TRUE){

        if (l == 1){
          beta_init <- complex((p-1), real = 0, imaginary = 0)
        }else{
          beta_init <- arr_beta[(l-1), ]
        }

        if (type == "ols"){

          arr_beta[l, ] <- classo_fixed(X = Xk, Y = yk,
                                        b.init = beta_init, lambda = Lambdas[l],
                                        weights = W[-c(k),k])
        }else if (type == "wls"){

          arr_beta[l, ] <- classo_cov_fixed(X = Xk, Y = yk,
                                            b.init = beta_init, lambda = Lambdas[l],
                                            weights = W[-c(k),k])
        }

      }else{

        if (type == "ols"){
          arr_beta[l, ] <- classo_fixed(X = Xk, Y = yk, lambda = Lambdas[l],
                                        scale, weights = W[-c(k),k])
        }else if (type == "wls"){
          arr_beta[l, ] <- classo_cov_fixed(X = Xk, Y = yk, lambda = Lambdas[l],
                                            scale, weights = W[-c(k),k])
        }

      }

      mse[l,k] <- sum(Mod(yk - Xk %*% arr_beta[l,])^2)/N
      bic_beta[l,k] <- log(mse[l,k]) + log(N) * sum(arr_beta[l, ]!=0)/N

    }

    idx_hat <- which.min(na.omit(bic_beta[,k]))
    b_hat[k,] <- arr_beta[idx_hat, ]
    # vec[-k] <- round(abs(b_hat)/sum(abs(b_hat)), 5)
    G_container[k,-c(k)] <- b_hat[k,]
  }

  Theta <- array(NA, dim=c(p,p))
  for(i1 in 1:p){
    for(i2 in 1:p){
      if (i1 <= i2){
        if (method == "or"){
          Theta[i1,i2] <- max(abs(G_container[i1,i2]),
                              abs(G_container[i2,i1]))
        }else{
          Theta[i1,i2] <- min(abs(G_container[i1,i2]),
                              abs(G_container[i2,i1]))
        }
      }else{
        Theta[i1,i2] <- Theta[i2,i1]
      }
    }
  }

  output <- list(Theta=Theta,MSE=mse,BIC=bic_beta,Lambda=Lambdas)
  return(output)
}
