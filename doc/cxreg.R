## ----setup-pdf, include=FALSE-------------------------------------------------
knitr::opts_chunk$set(
  fig.width  = 5,
  fig.height = 4,
  dpi        = 96      
)

## ----eval=FALSE, message=FALSE------------------------------------------------
# library(devtools)
# devtools::install_github("yk748/cxreg")

## -----------------------------------------------------------------------------
library(cxreg)

## -----------------------------------------------------------------------------
data(classo_example)
x <- classo_example$x
y <- classo_example$y

## -----------------------------------------------------------------------------
fit <- classo(x,y)

## -----------------------------------------------------------------------------
plot(fit, xvar="lambda", label=TRUE)

## -----------------------------------------------------------------------------
plot(fit, xvar="norm", label=TRUE)

## -----------------------------------------------------------------------------
plot(fit, xvar="dev", label=TRUE)

## -----------------------------------------------------------------------------
any(fit$lambda == 0.1)
coef(fit, s=0.1, exact=FALSE)

## -----------------------------------------------------------------------------
coef(fit, s=0.1, exact=TRUE, x=x, y=y)

## -----------------------------------------------------------------------------
set.seed(29)
nx <- array(rnorm(5*20), c(5,20)) + (1+1i) * array(rnorm(5*20), c(5,20))
for (j in 1:20) {
  nx[,j] <- nx[,j] / sqrt(mean(Mod(nx[,j])^2))
}
predict(fit, newx = nx, s = c(0.1, 0.05), type="response")

## -----------------------------------------------------------------------------
predict(fit, newx = nx, s = c(0.1, 0.05), type="coefficient")

## -----------------------------------------------------------------------------
predict(fit, newx = nx, s = c(0.1, 0.05), type="nonzero")

## -----------------------------------------------------------------------------
cvfit <- cv.classo(x,y,trace.it = 1)

## -----------------------------------------------------------------------------
print(cvfit)

## -----------------------------------------------------------------------------
plot(cvfit)

## -----------------------------------------------------------------------------
cvfit$lambda.min

## -----------------------------------------------------------------------------
coef(cvfit, s = "lambda.min")

## -----------------------------------------------------------------------------
predict(cvfit, newx = x[1:5,], s = "lambda.min")

## -----------------------------------------------------------------------------
data(cglasso_example)
f_hat <- cglasso_example$f_hat
n <- cglasso_example$n

## -----------------------------------------------------------------------------
fit_cglasso_I <- cglasso(S=f_hat,type="I",nobs=n)

## -----------------------------------------------------------------------------
fit_cglasso_II <- cglasso(S=f_hat,type="II",nobs=n, nlambda=30, stop_criterion = "AIC")

## -----------------------------------------------------------------------------
fit_cglasso_another <- cglasso(S=f_hat,type="II",nobs=n, stopping_rule = FALSE)
fit_cglasso_another$lambda_grid

## -----------------------------------------------------------------------------
plot(fit_cglasso_I$Theta_list,index=fit_cglasso_I$min_index,type="mod",label=TRUE)
plot(fit_cglasso_I$Theta_list,index=fit_cglasso_I$min_index,type="both",label=TRUE)
plot(fit_cglasso_II$Theta_list,index=fit_cglasso_II$min_index,type="real",label=FALSE)
plot(fit_cglasso_II$Theta_list,index=fit_cglasso_II$min_index,type="imaginary",label=FALSE)

## ----compact-pdf, include=FALSE, eval=TRUE------------------------------------
tryCatch({
  Sys.setenv(R_GSCMD = "C:/Program Files/gs/gs9.56.1/bin/gswin64c.exe")
  f <- tools::file_path_as_absolute("cxreg.pdf")
  if (file.exists(f)) {
    tools::compactPDF(paths = f, gs_quality = "ebook")
  }
}, error = function(e) invisible(NULL))

