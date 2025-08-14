---
title: 'cxreg: An R package for Complex-valued lasso and graphical lasso'
tags:
  - Spectral analysis
  - Regularization path
  - Coordinate descent 
  - lasso
  - Graphical model
authors:
  - name: Younghoon Kim
    orcid: 0009-0007-0117-5530
    equal-contrib: true
    affiliation: 1
  - name: Navonil Deb
    equal-contrib: true 
    affiliation: 1
  - name: Sumanta Basu
    corresponding: true 
    affiliation: 1
affiliations:
 - name: Department of Statistics and Data Science, Cornell University, United States
   index: 1
date: 2025-08-14
bibliography: paper.bib
---

# Summary

*Pathwise coordinate descent* [@friedman2007pathwise] is a popular iterative optimization method for solving penalized linear regression [*lasso*; @tibshirani1996regression] and penalized Gaussian log-likelihood [*graphical lasso*; e.g., @banerjee2008model] problems. Given a specified or derived sequence of tuning parameters $\lambda$ in the penalty term, the algorithm updates one coordinate at a time by isolating the target variable using partial residuals.

The coordinate descent algorithm for complex-valued lasso and Gaussian graphical lasso, along with the corresponding \texttt{cxreg} R package, leverages the standard pathwise coordinate descent method by demonstrating its applicability to complex-valued settings [@deb2024regularized]. Our algorithm is defined on the realifying complex numbers via a ring isomorphism [e.g., @herstein1991topics]. Specifically, for a complex number $z\in\mathbb{C}$, we define the map:

$$
\varphi(z) = \begin{pmatrix}
\mathbf{Re}(z) & -\mathbf{Im}(z) \\
\mathbf{Im}(z) & \mathbf{Re}(z)
\end{pmatrix}.
$$
We show that the $\varphi$ is a field isomorphism between $\mathbb{C}$ and the set

$$
\mathcal{M}^{2 \times 2} :=\bigg\{\begin{pmatrix} 
a & - b \\ b & a
\end{pmatrix} : a,b\in\mathbb{R}\bigg\}.
$$

Under this mapping, there exists connections between algebraic operations in the real and complex domains. We extend these operations to the $p$-dimensional case for $p>1$. Further details are provided in Section 3.1 of @deb2024regularized.

We use this mapping to transform linear regression problems in the complex-valued setting into equivalent problems in the real-valued setting. The same mapping can then be applied to the Gaussian log-likelihood for complex-valued variables. For the $\ell_1$  penalty term used in lasso and graphical lasso, we exploit the fact that the modulus of each complex-valued coordinate is equivalent to the entrywise $\ell_2$ norm of its real and imaginary parts. This implies that the lasso formulation for $p$-dimensional complex-valued variables is equivalent to a $2p$-dimensional group lasso [@yuan2006model], where each group consists of the real and imaginary components of a single complex variable ($p$ groups in total).


Finally, we implement the computational shortcuts, including *warm start* with the solutions for the sequence of the tuning parameter $\lambda$ in lasso and graphical lasso, and *active set selection* to reduce the number of predictors that need to be updated at each $\lambda$ [e.g., Chapter 5.4 in @hastie2015statistical]. In the R package, the user interface, such as the names of inputs for functions, is borrowed from the well-known pathwise coordinate descent R package \texttt{glmnet} [@friedman2010regularization]. A detailed vignette is available [Here](https://github.com/yk748/cxreg/blob/main/vignette/cxreg.pdf).


# Statement of need

Solving complex-valued penalized Gaussian likelihood (as well as penalized linear regression as a part of the problem) is critical in high-dimensional time series analysis, as examining partial spectral coherence in the frequency domain is analogous to examining partial correlations in Gaussian graphical models [e.g., @priestley1988spectral]. With the development of the local Whittle likelihood approximation [@whittle1951hypothesis], the computation of the inverse spectral density matrix, called the *spectral precision matrix*, enables the investigation of associations among variables in high-dimensional time series data.

These associations, captured by the spectral precision matrix, are evaluated at fixed frequencies. However, due to the inconsistency of the spectral density matrix estimator at a single frequency [e.g., Chapter 10 in @brockwell1991time], it is necessary to compute an averaged smoothed periodogram, obtained by aggregating periodogram matrices across neighboring frequencies, as the sample analog of the spectral density matrix. As a consequence, spectral precision matrices must be computed at multiple neighboring frequencies around each fixed frequency of interest. This is important in statistical inference of high-dimensional spectral precision matrix [@krampe2025frequency], requiring repetitive computations over high-dimensional objects. In this context, the development of efficient optimization algorithms for solving complex-valued penalized Gaussian likelihood is essential.


# State of the field

Although coordinate descent is known to be computationally efficient and has been widely used to solve the lasso [@friedman2007pathwise] and graphical Lasso [@friedman2008sparse], particularly under sparsity assumptions, it has not been extended to complex-valued settings. This limitation arises primarily because updating complex-valued partial residuals within the coordinate descent framework had not been explored.

The most well-known approach to this problem is based on the algorithm proposed in @fiecas2019spectral, which leverages *clime*  [@cai2011constrained], formulated via linear programming. An alternative is the use of the alternating direction method of multipliers (*ADMM*), as described in @baek2023local. However, neither of these methods explicitly exploits the sparsity nature of the model, even though the problem is defined in high-dimensional regimes. Our proposed algorithm is designed to be highly efficient in sparse settings, which has been compared to the benchmarks in terms of both estimation accuracy and computational efficiency. The details can be found in @deb2024regularized.


# Acknowledgements

ND and SB acknowledge partial support from NIH award R21NS120227. In addition, SB acknowledges partial support from NSF awards DMS-1812128, DMS-2210675, DMS-2239102 and NIH award R01GM135926.

# References

