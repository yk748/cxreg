
# Complex-valued Lasso

We provide efficient procedures for fitting the entire complex-valued
lasso, called classo, for linear regression. The idea is to bring the
exact coordinate descent algorithm ([2010](#ref-glmnet)) into complex
analog, by introducing a unique isomorphism. Details may be found in
Deb, Kuceyeski, and Basu ([2013](#ref-classo)).

Currently, the beta version (0.1.1) is released, and the version
provides the basic functions, including fit, cross-validation, and other
functions for illustration. We will update the error control and other
functions for users’ convenience accordingly. Please email Younghoon Kim
<yk748@cornell.edu> if any bugs/errors have been discovered.

## References

<div id="refs" class="references">

<div id="ref-classo">

Navonil Deb, Amy Kuceyeski, and Sumanta Basu. 2024. “Regularized
Estimation of Sparse Spectral Precision Matrices.” *arXiv preprint
arXiv:2401.11128.*. <https://arxiv.org/abs/2401.11128>.

</div>

<div id="refs" class="references">

<div id="ref-glmnet">

Friedman, Jerome, Trevor Hastie, and Robert Tibshirani. 2010.
“Regularization Paths for Generalized Linear Models via Coordinate
Descent.” *Journal of Statistical Software, Articles* 33 (1): 1–22.
<https://doi.org/10.18637/jss.v033.i01>.

</div>

</div>

</div>
