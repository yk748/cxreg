
# Complex-valued lasso (CLASSO) and complex-valued graphical lasso (CGLASSO)

We provide an efficient estimation procedure for fitting a penalized
complex-valued lasso (CLASSO) and a complex-valued graphical lasso
(CGLASSO) problem. We implement a pathwise coordinate descent algorithm
(Friedman et al. 2007; Friedman, Hastie, and Tibshirani 2008; Friedman,
Hastie, and Tibshirani 2010) into a complex analog, by introducing an
isomorphism between complex numbers and a set of $2\times 2$ orthogonal
matrices. Such an isomorphism enables leveraging the existing algorithms
in `glmnet` for complex settings. Details can be found in Deb,
Kuceyeski, and Basu (2024).

## `cxreg` installation for macOS

Following are the instructions to install the GCC compiler that includes
`gfortran` for using in `R`.

### GCC and `gfortran` installation steps for `R`

Open terminal and run:

``` bash
brew install gcc
```

Check if the `/.R` directory already exists by running:

``` bash
ls -ld ~/.R
```

If the `/.R` does not exist run the first line. Open the `/.R/Makevars`
file by running the second line.

``` bash
mkdir ~/.R
vim ~/.R/Makevars
```

Add the following lines in `/.R/Makevars` and save the file.

``` bash
FC = /opt/homebrew/Cellar/gcc/15.0.1/bin/gfortran
F77 = /opt/homebrew/Cellar/gcc/15.0.1/bin/gfortran
FLIBS = -L/opt/homebrew/Cellar/gcc/15.0.1/lib/gcc/15
```

One needs to change the GCC version `15.0.1` to the according version
obtained by `gfortran --version`.

### Package installation steps

Install the package from GitHub using `devtools` as follows:

``` r
library(devtools)
devtools::install_github("yk748/cxreg")
```

A short example is presented to demonstrate the working of the package.
Details can be found in `cxreg/vignettes/cxreg.Rmd`.

``` r
library(cxreg)
data(classo_example)
x <- classo_example$x
y <- classo_example$y
cvfit <- cv.classo(x,y,trace.it = 1)
print(coef(cvfit, s = "lambda.min"))
predict(cvfit, newx = x[1:5,], s = "lambda.min")
```

The first stable version (1.0.0) is released, which provides the fitting
functions of CLASSO and CGLASSO along with their auxiliary functions
such as printing paths of coefficients, cross-validation, and generating
plots (regularization paths and heatmap). We will keep updating and
maintaining the package to eventually incorporate a systematic error
control for users’ convenience accordingly. Please email Younghoon Kim
<yk748@cornell.edu> if any bugs/errors have been discovered. All
remaining errors are our own.

## References

Navonil Deb, Amy Kuceyeski, and Sumanta Basu. 2024. “Regularized
Estimation of Sparse Spectral Precision Matrices.” *arXiv preprint
arXiv:2401.11128.*. <https://arxiv.org/abs/2401.11128>.

Jerome Friedman, Trevor Hastie, Holger Hofling, and Robert Tibshirani.
2007. “Pathwise Coordinate Optimization.” *The Annals of Applied
Statistics* 1(2): 302-332.
<https://web.archive.org/web/20170301123147id_/http://gautampendse.com/software/lasso/webpage/Friedman2007.pdf>

Friedman, Jerome, Trevor Hastie, and Robert Tibshirani. 2008. “Sparse
inverse covariance estimation with the graphical lasso.” *Biostatistics*
9(3): 432-441.
<https://www.asc.ohio-state.edu/statistics/statgen/joul_aut2015/2008-Friedman-Hastie-Tibshirani.pdf>

Jerome Friedman, Trevor Hastie, and Robert Tibshirani. 2010.
“Regularization Paths for Generalized Linear Models via Coordinate
Descent.” *Journal of Statistical Software, Articles* 33 (1): 1–22.
<https://doi.org/10.18637/jss.v033.i01>.
