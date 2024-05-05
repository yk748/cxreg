rm(list=ls())
library(pROC)
library(plot3D)
library(gdata)
library(clime)
library(mvtnorm)
library(gplots)
library(superheat)

source("Classo_fixed_YK.R")
source("Classo_cov_fixed_YK.R")
source("Complex_operators_YK.R")
source("library_misc_YK.R")
source("Nodewise_YK.R")
# ------------------------------------------------ #
# Data generation
# ------------------------------------------------ #
# Size of problem
p <- 30
n <- 500
burn <- 2*n

# ------------------------------------------------ #
# parameter settings
# DGP 1
a <- 0.7; b <- 0.3
C <- array(0,dim=c(p,p))
for (i in 1:p){
  C[i,i] <- a
  
  if (i < p){
    C[i,(i+1)] <- b
  }
  if (i > 1){
    C[i,(i-1)] <- b
  }
}
Sigma <- solve(C)
f_dgp1 <- array(NA,dim=c(p,p,n))
for (j in c(0:(n-1))){
  f[,,(j+1)] <- Sigma /(2*pi)
}
# 
# 
# # DGP 2
# a1 <- 0.7; b1 <- 0.3
# a2 <- 0.2; b2 <- 0
# C1 <- array(0,dim=c(p,p))
# C2 <- array(0,dim=c(p,p))
# for (i in 1:p){
#   C1[i,i] <- a1
#   C2[i,i] <- a2
#   
#   if (i < p){
#     C1[i,(i+1)] <- b1
#     C2[i,(i+1)] <- b2
#   }
#   if (i > 1){
#     C1[i,(i-1)] <- b1
#     C2[i,(i-1)] <- b2
#   }
# }
# Sigma1 <- solve(C1)
# Sigma2 <- solve(C2)

# ------------------------------------------------ #
# Compute spectral densities and PSC
modulo <- seq(0,2*pi*(n-1)/n,by=2*pi/n)
j <- 1
Theta_j <- solve(f[,,j])
D_f_j <- diag(diag(Theta_j),p)
PSC_j <- solve(sqrt(D_f_j)) %*% Theta_j %*% solve(sqrt(D_f_j)) 

# ------------------------------------------------ #
# Sampling
Xt <- rmvnorm(n = n, mean = rep(0, p), sigma = Sigma)
# Xt <- rbind(rmvnorm(n = n/2, mean = rep(0, p), sigma = Sigma1),
#             rmvnorm(n = n/2, mean = rep(0, p), sigma = Sigma2))
# plot.ts(Xt[,1])

# ------------------------------------------------ #
# Data generation
# ------------------------------------------------ #
Theta_j_ols <- nodewise_YK(Xt, j, method="or", type = "ols", pathwise=FALSE)
Theta_j_ols_path <- nodewise_YK(Xt, j, method="or", type = "ols", pathwise=TRUE)
Theta_j_wls <- nodewise_YK(Xt, j, method="or", type = "wls", pathwise=FALSE)
Theta_j_wls_path <- nodewise_YK(Xt, j, method="or", type = "wls", pathwise=TRUE)

# # Under construction:
# Theta_hat_0 <- nodewise_YK(Xt, j, method="or", type = "ols", 
#                            scale = FALSE,  pathwise=FALSE)$Theta
# weight_mat <- diag(Theta_hat_0) %*% t(diag(Theta_hat_0))
# Theta_j_ols_adapt <- nodewise_YK(Xt, j, method="or", type = "ols", 
#                                  scale = FALSE, weights=weight_mat, 
#                                  pathwise=FALSE)

# ------------------------------------------------ #
# Comparison
# ------------------------------------------------ #
superheat(PSC_j,heat.col.scheme = "red")
superheat(Theta_j_ols$Theta,heat.col.scheme = "red")
superheat(Theta_j_ols_path$Theta,heat.col.scheme = "red")
superheat(Theta_j_wls$Theta,heat.col.scheme = "red")
superheat(Theta_j_wls_path$Theta,heat.col.scheme = "red")
