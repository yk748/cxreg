#' dft.X
#'
#' discrete Fourier transformation
#' @param X input1
#' @param j input2
#' @param m input3
#' @return output
#' @examples
#' example <- dft.X(X, j, m)
#'
#' @export
dft.X <- function(X, j, m){
  n <- nrow(X)
  p <- ncol(X)

  dft <- mvfft(X)/sqrt(2*pi*n)

  vj <- j + c(-m:m)
  ind <- fixm(vj, n)
  Z <- dft[ind, ]

  return(Z)
}

#' fixm
#'
#' window
#' @param v input1
#' @param n input2
#' @return output
#' @examples
#' example <- fixm (v, n)
#'
#' @export
fixm <- function(v, n){

  v[v<1] = v[v<1]+n
  v[v>n] = v[v>n]-n

  return(v)
}

# -------------------------------------------- #
