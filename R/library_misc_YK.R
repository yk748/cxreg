# -------------------------------------------- #
dft.X <- function(X, j, m){
  n <- nrow(X)
  p <- ncol(X)
  
  dft <- mvfft(X)/sqrt(2*pi*n)
  
  vj <- j + c(-m:m)
  ind <- fixm(vj, n)
  Z <- dft[ind, ]
  
  return(Z)
}

# -------------------------------------------- #
fixm <- function(v, n){
  
  v[v<1] = v[v<1]+n
  v[v>n] = v[v>n]-n
  
  return(v)
}

# -------------------------------------------- #