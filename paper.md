# Summary

*Pathwise coordinate descent* (Friedman et al. 2007) is a popular
iterative optimization method for solving penalized linear regression
\[*lasso*; Tibshirani (1996)\] and penalized Gaussian log-likelihood
\[*graphical lasso*; e.g., Banerjee, El Ghaoui, and d’Aspremont (2008)\]
problems. Given a specified or derived sequence of tuning parameters
$\lambda$ in the penalty term, the algorithm updates one coordinate at a
time by isolating the target variable using partial residuals.

The coordinate descent algorithm for complex-valued lasso and Gaussian
graphical lasso, along with the corresponding R package, leverages the
standard pathwise coordinate descent method by demonstrating its
applicability to complex-valued settings (Deb, Kuceyeski, and Basu
2024). Our algorithm is defined on the realifying complex numbers via a
ring isomorphism (e.g., Herstein 1991; Dummit, Foote, et al. 2004).
Specifically, for a complex number $z\in\mathbb{C}$, we define the map:

$$
\varphi(z) = \begin{pmatrix}
\mathbf{Re}(z) & -\mathbf{Im}(z) \\
\mathbf{Im}(z) & \mathbf{Re}(z)
\end{pmatrix}.
$$ We show that the $\varphi$ is a field isomorphism between
$\mathbb{C}$ and the set

$$
\mathcal{M}^{2 \times 2} :=\left\{\begin{pmatrix} 
a & - b \\ b & a
\end{pmatrix} : a,b\in\mathbb{R}\right\}.
$$

Under this mapping, there exists connections between algebraic
operations in the real and complex domains. We extend these operations
to the $p$-dimensional case for $p>1$. Further details are provided in
Section 3.1 of Deb, Kuceyeski, and Basu (2024).

We use this mapping to transform linear regression problems in the
complex-valued setting into equivalent problems in the real-valued
setting. The same mapping can then be applied to the Gaussian
log-likelihood for complex-valued variables. For the $\ell_1$ penalty
term used in lasso and graphical lasso, we exploit the fact that the
modulus of each complex-valued coordinate is equivalent to the entrywise
$\ell_2$ norm of its real and imaginary parts. This implies that the
lasso formulation for $p$-dimensional complex-valued variables is
equivalent to a $2p$-dimensional group lasso (Yuan and Lin 2006), where
each group consists of the real and imaginary components of a single
complex variable ($p$ groups in total).

Finally, we implement the computational shortcuts, including *warm
start* with the solutions for the sequence of the tuning parameter
$\lambda$ in lasso and graphical lasso, and *active set selection* to
reduce the number of predictors that need to be updated at each
$\lambda$ (e.g., Chapter 5.4 in Hastie, Tibshirani, and Wainwright
2015). In the R package, the user interface, such as the names of inputs
for functions, is borrowed from the well-known pathwise coordinate
descent R package (Friedman, Hastie, and Tibshirani 2010). A detailed
vignette is available
[Here](https://github.com/yk748/cxreg/blob/main/vignettes/cxreg.pdf).

# Statement of need

Solving complex-valued penalized Gaussian likelihood (as well as
penalized linear regression as a part of the problem) is critical in
high-dimensional time series analysis, as examining partial spectral
coherence in the frequency domain is analogous to examining partial
correlations in Gaussian graphical models (e.g., Priestley 1988;
Dahlhaus 2000). With the development of the local Whittle likelihood
approximation (Whittle 1951), the computation of the inverse spectral
density matrix, called the *spectral precision matrix*, enables the
investigation of associations among variables in high-dimensional time
series data.

These associations, captured by the spectral precision matrix, are
evaluated at fixed frequencies. However, due to the inconsistency of the
spectral density matrix estimator at a single frequency (e.g., Chapter
10 in Brockwell and Davis 1991), it is necessary to compute an averaged
smoothed periodogram, obtained by aggregating periodogram matrices
across neighboring frequencies, as the sample analog of the spectral
density matrix. As a consequence, spectral precision matrices must be
computed at multiple neighboring frequencies around each fixed frequency
of interest. This is important in statistical inference of
high-dimensional spectral precision matrix (Krampe and Paparoditis
2025), requiring repetitive computations over high-dimensional objects.
In this context, the development of efficient optimization algorithms
for solving complex-valued penalized Gaussian likelihood is essential.

# State of the field

Although coordinate descent is known to be computationally efficient and
has been widely used to solve the lasso (Friedman et al. 2007) and
graphical Lasso (Friedman, Hastie, and Tibshirani 2008), particularly
under sparsity assumptions, it has not been extended to complex-valued
settings. This limitation arises primarily because updating
complex-valued partial residuals within the coordinate descent framework
had not been explored.

The most well-known approach to this problem is based on the algorithm
proposed in Fiecas et al. (2019), which leverages *clime* (Cai, Liu, and
Luo 2011), formulated via linear programming. An alternative is the use
of the alternating direction method of multipliers (*ADMM*), as
described in Baek, Düker, and Pipiras (2023). However, neither of these
methods explicitly exploits the sparsity nature of the model, even
though the problem is defined in high-dimensional regimes. Our proposed
algorithm is designed to be highly efficient in sparse settings, which
has been compared to the benchmarks in terms of both estimation accuracy
and computational efficiency. The details can be found in Deb,
Kuceyeski, and Basu (2024).

# Acknowledgements

ND and SB acknowledge partial support from NIH award R21NS120227. In
addition, SB acknowledges partial support from NSF awards DMS-1812128,
DMS-2210675, DMS-2239102 and NIH award R01GM135926.

# References

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-baek2023local" class="csl-entry">

Baek, Changryong, Marie-Christine Düker, and Vladas Pipiras. 2023.
“Local Whittle Estimation of High-Dimensional Long-Run Variance and
Precision Matrices.” *The Annals of Statistics* 51 (6): 2386–2414.

</div>

<div id="ref-banerjee2008model" class="csl-entry">

Banerjee, Onureena, Laurent El Ghaoui, and Alexandre d’Aspremont. 2008.
“Model Selection Through Sparse Maximum Likelihood Estimation for
Multivariate Gaussian or Binary Data.” *The Journal of Machine Learning
Research* 9: 485–516.

</div>

<div id="ref-brockwell1991time" class="csl-entry">

Brockwell, Peter J, and Richard A Davis. 1991. *Time Series: Theory and
Methods*. Springer science & business media.

</div>

<div id="ref-cai2011constrained" class="csl-entry">

Cai, Tony, Weidong Liu, and Xi Luo. 2011. “A Constrained $\ell_1$
Minimization Approach to Sparse Precision Matrix Estimation.” *Journal
of the American Statistical Association* 106 (494): 594–607.

</div>

<div id="ref-dahlhaus2000graphical" class="csl-entry">

Dahlhaus, Rainer. 2000. “Graphical Interaction Models for Multivariate
Time Series1.” *Metrika* 51 (2): 157–72.

</div>

<div id="ref-deb2024regularized" class="csl-entry">

Deb, Navonil, Amy Kuceyeski, and Sumanta Basu. 2024. “Regularized
Estimation of Sparse Spectral Precision Matrices.” *arXiv Preprint
arXiv:2401.11128*.

</div>

<div id="ref-dummit2004abstract" class="csl-entry">

Dummit, David Steven, Richard M Foote, et al. 2004. *Abstract Algebra*.
Vol. 3. Wiley Hoboken.

</div>

<div id="ref-fiecas2019spectral" class="csl-entry">

Fiecas, Mark, Chenlei Leng, Weidong Liu, and Yi Yu. 2019. “Spectral
Analysis of High-Dimensional Time Series.” *Electronic Journal of
Statistics* 13 (2): 4079–4101.

</div>

<div id="ref-friedman2007pathwise" class="csl-entry">

Friedman, Jerome, Trevor Hastie, Holger Höfling, and Robert Tibshirani.
2007. “Pathwise Coordinate Optimization.” *The Annals of Applied
Statistics*, 302–32.

</div>

<div id="ref-friedman2010regularization" class="csl-entry">

Friedman, Jerome, Trevor Hastie, and Rob Tibshirani. 2010.
“Regularization Paths for Generalized Linear Models via Coordinate
Descent.” *Journal of Statistical Software* 33: 1–22.

</div>

<div id="ref-friedman2008sparse" class="csl-entry">

Friedman, Jerome, Trevor Hastie, and Robert Tibshirani. 2008. “Sparse
Inverse Covariance Estimation with the Graphical Lasso.” *Biostatistics*
9 (3): 432–41.

</div>

<div id="ref-hastie2015statistical" class="csl-entry">

Hastie, Trevor, Robert Tibshirani, and Martin Wainwright. 2015.
“Statistical Learning with Sparsity.” *Monographs on Statistics and
Applied Probability* 143 (143): 8.

</div>

<div id="ref-herstein1991topics" class="csl-entry">

Herstein, Israel N. 1991. *Topics in Agebra*. John Wiley & Sons.

</div>

<div id="ref-krampe2025frequency" class="csl-entry">

Krampe, Jonas, and Efstathios Paparoditis. 2025. “Frequency Domain
Statistical Inference for High-Dimensional Time Series.” *Journal of the
American Statistical Association*, 1–13.

</div>

<div id="ref-priestley1988spectral" class="csl-entry">

Priestley, Maurice Bertram. 1988. *The Spectral Analysis of Time
Series*. Oxford University Press.

</div>

<div id="ref-tibshirani1996regression" class="csl-entry">

Tibshirani, Robert. 1996. “Regression Shrinkage and Selection via the
Lasso.” *Journal of the Royal Statistical Society Series B: Statistical
Methodology* 58 (1): 267–88.

</div>

<div id="ref-whittle1951hypothesis" class="csl-entry">

Whittle, Peter. 1951. “Hypothesis Testing in Time Series Analysis.” *(No
Title)*.

</div>

<div id="ref-yuan2006model" class="csl-entry">

Yuan, Ming, and Yi Lin. 2006. “Model Selection and Estimation in
Regression with Grouped Variables.” *Journal of the Royal Statistical
Society Series B: Statistical Methodology* 68 (1): 49–67.

</div>

</div>
