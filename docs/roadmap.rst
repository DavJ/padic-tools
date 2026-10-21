Scope, release readiness and roadmap
====================================

Version 0.1 is a working research library with explicit mathematical domains.
It is not a complete replacement for Sage, PARI/GP or a general computer
algebra system. A useful notion of completeness is a documented API with
verified behavior in its declared scope, followed by staged extensions.

Implemented
-----------

* Exact rational conversion and finite precision Qp ball arithmetic, including
  nonintegral values, valuation and precision checks.
* Local exp/log/sin/cos, principal-unit powers, odd/dyadic square roots and
  simple-root Hensel lifting.
* Multiplicative theta3 and theta4 constants with tail bounds.
* Exact Bernoulli/classical-zeta special values and numerical even KL branches
  on Zp away from the point s=1, with a documented branch convention.
* Exact rational polynomials and formal power series; optional SymPy function
  nodes and numerical evaluation with domain checks.
* Rational/Gaussian diagonal bridges, selected split places, quaternion and
  biquaternion scalar base changes, and an explicit even-Clifford isomorphism.
* A source-linked HTML guide and executable examples.

Next mathematical features
--------------------------

1. Broader independent Sage/PARI oracle coverage and performance baselines,
   particularly around KL branches and near poles.
2. Ramified/unramified extension fields with uniformizers, residue fields and
   extension-aware valuation/precision propagation. The current quadratic
   coefficient algebra remains useful but is not that full engine.
3. Singular Hensel lifting, general polynomial factorization, matrix solvers
   with condition-aware precision and optional accelerated arithmetic backends.
4. General Dirichlet characters, the removable s=1 limits of nonprincipal
   branches, p-adic gamma/polylogarithms and larger-prime zeta algorithms.
5. Theta1/theta2 with explicit fractional-power choices, lattice theta series
   and the geometric data required for rigid-analytic theta comparisons.
6. General number-field places and invariants that have proved comparison
   theorems across completions. Any proposed complex-zero bridge needs a
   separate theorem and must not be inferred from agreement at special values.

PyPI release preparation
------------------------

``pyproject.toml`` defines the distribution, import package, optional extras,
metadata and build backend. Both wheel and source distribution are generated
with ``python -m build``. The release workflow is manual and defaults to
TestPyPI; no upload is triggered by an ordinary push or pull request.

Before a first public release the owner must select a license and add its
text/metadata, settle the distribution/import name, configure a PyPI account
and Trusted Publisher, and approve the actual upload. A missing license is a
project release-policy gate in the supplied workflow, not a claim that PyPI
technically refuses all unlicensed distributions.

No new license has been assigned to the owner's existing work. The manual
contains third-party references rather than copied implementations. Existing
package-name availability is only a point-in-time observation; a successful
release, not a failed lookup, establishes that a name can be registered.
