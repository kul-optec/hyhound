[![arXiv Preprint](https://img.shields.io/badge/arXiv-Preprint-b31b1b)](https://arxiv.org/abs/2503.15372v1)
[![CI: Linux](https://github.com/kul-optec/hyhound/actions/workflows/linux.yml/badge.svg)](https://github.com/kul-optec/hyhound/actions/workflows/linux.yml)
[![PyPI Downloads](https://img.shields.io/pypi/dm/hyhound?label=PyPI&logo=python)](https://pypi.org/project/hyhound)

# hyhound

**Hy**perbolic **Ho**useholder transformations for **U**p- a**n**d **D**owndating Cholesky factorizations.

## Purpose

Given a Cholesky factor $L$ of a dense matrix $H$, the `hyhound::update_cholesky`
function computes the Cholesky factor $\tilde L$ of the matrix
$\tilde H = \tilde L \tilde L^\top = H + A \Sigma A^\top$ (where
$H,\tilde H\in\mathbb{R}^{n\times n}$ with $H \succ 0$ and $\tilde H \succ 0$,
$A \in \mathbb{R}^{n\times m}$, $\Sigma \in \mathbb{R}^{m\times m}$ diagonal,
and $L, \tilde L\in\mathbb{R}^{n\times n}$ lower triangular).

Computing $\tilde L$ in this way is done in $mn^2 + \mathcal{O}(n^2 + mn)$
operations rather than the $\frac16 n^3 + \frac12 mn^2 + \mathcal{O}(n^2 + mn)$
operations required for the explicit evaluation and factorization of $\tilde H$.
When $m \ll n$, this results in a considerable speedup over full factorization,
enabling efficient low-rank updates of Cholesky factorizations, for use in e.g.
iterative algorithms for numerical optimization.

Additionally, hyhound includes efficient routines for updating factorizations
of the Riccati recursion for optimal control problems.

## Preprint

The paper describing the algorithms in this repository can be found on arXiv: <https://arxiv.org/abs/2503.15372v1>

```bibtex
@misc{pas_blocked_2025,
	title = {Blocked {Cholesky} factorization updates of the {Riccati} recursion using hyperbolic {Householder} transformations},
	url = {http://arxiv.org/abs/2503.15372},
	doi = {10.48550/arXiv.2503.15372},
	publisher = {arXiv},
	author = {Pas, Pieter and Patrinos, Panagiotis},
	month = mar,
	year = {2025},
	note = {Accepted for publication in the Proceedings of CDC 2025}
}
```
