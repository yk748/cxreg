---
title: 'cxreg: An R package for complex-valued lasso and graphical lasso'
tags:
  - Spectral analysis
  - Regularization path
  - Coordinate descent 
  - Lasso
  - Graphical model
authors:
  - name: Younghoon Kim
    orcid: 0009-0007-0117-5530
    equal-contrib: true
    affiliation: 1
  - name: Navonil Deb
    equal-contrib: true 
    affiliation: 2
  - name: Sumanta Basu
    corresponding: true 
    affiliation: 2
affiliations:
 - name: Center for Data Science for Enterprise and Society, Cornell University, United States
   index: 1
 - name: Department of Statistics and Data Science, Cornell University, United States
   index: 2
date: 2026-04-20
citeproc: true
bibliography: paper.bib
---

# Summary

*Pathwise coordinate descent* [@friedman2007pathwise] is a widely used iterative optimization method for solving penalized linear regression (*lasso*; [@tibshirani1996regression]) and penalized Gaussian log-likelihood (*graphical lasso*; e.g., Banerjee, El Ghaoui, and d’Aspremont [@banerjee2008model]) problems. For a given sequence of tuning parameters $\lambda$ in the penalty term, the algorithm updates one coefficient at a time by computing partial residuals.

The coordinate descent algorithm for complex-valued lasso and Gaussian graphical lasso, implemented in the cxreg R package, extends the standard pathwise coordinate descent method to complex-valued data (Deb, Kuceyeski, and Basu [@deb2024regularized]). The algorithm operates by *realifying* complex numbers through a ring isomorphism [@herstein1991topics]. Specifically, for a complex number $z \in \mathbb{C}$, we define the map:

$$
\varphi(z) = 
\begin{pmatrix} 
\mathrm{Re}(z) & -\mathrm{Im}(z) \\ 
\mathrm{Im}(z) & \mathrm{Re}(z) 
\end{pmatrix}.
$$ 

It is known that $\varphi$ is a field isomorphism between
$\mathbb{C}$ and the set $\mathcal{M}^{2 \times 2}$ whose entries consists of $a,b \in \mathbb{R}$ to form 2 by 2 matrices [@artin2018algebra],

$$
\begin{pmatrix}
a & -b \\ 
b & a
\end{pmatrix}.
$$

Under this mapping, algebraic operations in the complex domain correspond directly to those in the real domain. We extend these operations to the $p$-dimensional case $(p>1)$; see Section 3.1 of Deb, Kuceyeski, and Basu [@deb2024regularized] for details.

Using this correspondence, we transform complex-valued linear regression problems into equivalent real-valued formulations. The same idea applies to the Gaussian log-likelihood for complex-valued variables. For the $\ell_{1}$ penalty used in lasso and graphical lasso, we note that the modulus of a complex coefficient equals the entrywise $\ell_{2}$ norm of its real and imaginary parts. Consequently, the lasso for $p$-dimensional complex variables is equivalent to a $2p$-dimensional group lasso (Yuan and Lin [@yuan2006model]), where each group contains the real and imaginary components of one complex variable.

We also incorporate standard computational enhancements, including warm starts, using solutions from previous $\lambda$ values, and active set selection to limit updates to relevant predictors (see Hastie, Tibshirani, and Wainwright [@hastie2015statistical], Chapter 5.4). The R package `glmnet` adopts a familiar interface similar to glmnet (Friedman, Hastie, and Tibshirani [@friedman2010regularization]), and a detailed vignette is available [Here](https://github.com/yk748/cxreg/blob/main/doc/cxreg.pdf).


# Statement of need

Solving complex-valued penalized Gaussian likelihood problems (and, as a component, penalized linear regression) is central to high-dimensional time series analysis. In the frequency domain, examining partial spectral coherence is analogous to studying partial correlations in Gaussian graphical models (e.g., Priestley [@priestley1988spectral]). With the development of the local Whittle likelihood approximation (Whittle [@whittle1951hypothesis]), the inverse spectral density matrix, known as the spectral precision matrix, provides a way to characterize associations among variables in high-dimensional time series data.

These associations, represented by the spectral precision matrix, are evaluated at fixed frequencies. However, because the spectral density matrix estimator is inconsistent at a single frequency (see Brockwell and Davis [@brockwell1991time], Chapter 10), it is common to use an averaged, smoothed periodogram, obtained by aggregating periodogram matrices across neighboring frequencies, as a consistent sample analog. Consequently, spectral precision matrices must be computed across multiple neighboring frequencies around each frequency of interest. Such analyses involve repeated high-dimensional optimization, as required for statistical inference on spectral precision matrices (Krampe and Paparoditis [@krampe2025frequency]). Therefore, developing efficient algorithms for complex-valued penalized Gaussian likelihood estimation is crucial for scalable and accurate frequency-domain analysis.


# State of the field

Although coordinate descent is computationally efficient and widely used to solve the lasso (Friedman et al. [@friedman2007pathwise]) and graphical lasso (Friedman, Hastie, and Tibshirani [@friedman2008sparse]), particularly under sparsity assumptions, its extension to complex-valued settings remains limited. This gap largely arises from the lack of methods for updating complex-valued partial residuals within the coordinate descent framework.

Several alternative approaches have been proposed for complex-valued graphical modeling, but many do not explicitly exploit sparsity within a coordinate descent framework, which can limit scalability in high-dimensional settings. Recent software, such as the R package tsglasso (Dallakyan et al. [@dallakyan2022time]) and the Python package FreDom (Dallakyan [@dallakyan2024learning]), provide related implementations based on closed-form updates derived via Wirtinger calculus. Another line of work employs the alternating direction method of multipliers (ADMM; Baek, Düker, and Pipiras [@baek2023local]).

In contrast, the algorithm implemented in this package is specifically designed to leverage sparsity within a coordinate descent framework. The method follows [@deb2024regularized], where comparisons with benchmark approaches, including the complex-valued extension of the nodewise regression CIPE algorithm (Fiecas et al. [@fiecas2019spectral]) and its pathwise variants, demonstrate improvements in both estimation accuracy and computational efficiency.


# Acknowledgements

ND and SB acknowledge partial support from NIH award R21NS120227. In addition, SB acknowledges partial support from NSF awards DMS-1812128, DMS-2210675, DMS-2239102 and NIH award R01GM135926. The authors thank the two referees for their helpful comments and suggestions, which improved the clarity and presentation of the manuscript.
