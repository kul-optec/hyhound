.. image:: https://img.shields.io/badge/arXiv-Preprint-b31b1b
   :target: https://arxiv.org/abs/2503.15372v1
   :alt: arXiv Preprint

.. image:: https://github.com/kul-optec/hyhound/actions/workflows/linux.yml/badge.svg
   :target: https://github.com/kul-optec/hyhound/actions/workflows/linux.yml
   :alt: CI: Linux

.. image:: https://img.shields.io/pypi/dm/hyhound?label=PyPI&logo=python
   :target: https://pypi.org/project/hyhound
   :alt: PyPI Downloads


hyhound
=======

**Hy**\perbolic **Ho**\useholder transformations for **U**\p- ‘**n**’ **D**\owndating Cholesky factorizations.


Purpose
-------

Given a Cholesky factor :math:`L` of a dense matrix :math:`H`, the
``hyhound::update_cholesky`` function computes the Cholesky factor
:math:`\tilde L` of the matrix

.. math::

   \tilde H = \tilde L \tilde L^\top = H + A \Sigma A^\top,

where :math:`H,\tilde H\in\mathbb{R}^{n\times n}` with :math:`H \succ 0`
and :math:`\tilde H \succ 0`, :math:`A \in \mathbb{R}^{n\times m}`,
:math:`\Sigma \in \mathbb{R}^{m\times m}` diagonal,
and :math:`L, \tilde L\in\mathbb{R}^{n\times n}` lower triangular.

Computing :math:`\tilde L` in this way is done in
:math:`mn^2 + \mathcal{O}(n^2 + mn)` operations rather than the
:math:`\tfrac16 n^3 + \tfrac12 mn^2 + \mathcal{O}(n^2 + mn)` operations
required for the explicit evaluation and factorization of :math:`\tilde H`.
When :math:`m \ll n`, this results in a considerable speedup over full
factorization, enabling efficient low-rank updates of Cholesky
factorizations, for use in e.g. iterative algorithms for numerical
optimization.

Additionally, hyhound includes efficient routines for updating
factorizations of the Riccati recursion for optimal control problems.


Preprint
--------

The paper describing the algorithms in this repository can be found on arXiv:  
`https://arxiv.org/abs/2503.15372v1 <https://arxiv.org/abs/2503.15372v1>`_

.. code-block:: bibtex

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
