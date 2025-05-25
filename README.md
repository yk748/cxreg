
# Complex-valued Lasso and graphical Lasso

We provide an efficient estimation procedure for fitting a penalized
complex-valued Lasso fitting (classo) and a complex-valued graphical
Lasso (cglasso) models. The idea is to bring the pathwise coordinate
descent algorithm
([2007](#ref-pathwise),[2008](#ref-glasso),[2010](#ref-glmnet)) into
complex analog, by introducing a unique isomorphism. Details can be
found in Deb, Kuceyeski, and Basu ([2024](#ref-classo)).

Currently, the beta version (0.8.0) is released, and the version
provides the basic fitting functions and their auxiliary functions,
including printing paths of coefficients, cross-validation, and
generating plots. We will keep updating and maintaining the package to
eventually incorporate a systematic error control for users’ convenience
accordingly. Please email Younghoon Kim <yk748@cornell.edu> if any
bugs/errors have been discovered.

## References

<div id="ref-classo" class="references">

Navonil Deb, Amy Kuceyeski, and Sumanta Basu. 2024. “Regularized
Estimation of Sparse Spectral Precision Matrices.” *arXiv preprint
arXiv:2401.11128.*. <https://arxiv.org/abs/2401.11128>.

<div id="ref-pathwise" class="references">

Jerome Friedman, Trevor Hastie, Holger Hofling, and Robert Tibshirani.
2007. “Pathwise Coordinate Optimization.” *The Annals of Applied
Statistics* 1(2): 302-332.
<https://web.archive.org/web/20170301123147id_/http://gautampendse.com/software/lasso/webpage/Friedman2007.pdf>

<div id="ref-glasso" class="references">

Friedman, Jerome, Trevor Hastie, and Robert Tibshirani. 2008. “Sparse
inverse covariance estimation with the graphical lasso.” *Biostatistics*
9(3): 432-441.
<https://www.asc.ohio-state.edu/statistics/statgen/joul_aut2015/2008-Friedman-Hastie-Tibshirani.pdf>

<div id="ref-glmnet" class="references">

Jerome Friedman, Trevor Hastie, and Robert Tibshirani. 2010.
“Regularization Paths for Generalized Linear Models via Coordinate
Descent.” *Journal of Statistical Software, Articles* 33 (1): 1–22.
<https://doi.org/10.18637/jss.v033.i01>.

