Architecture and contribution guide
===================================

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Module
     - Responsibility
   * - ``core``
     - Exact rational conversion, prime validation and immutable Qp balls
   * - ``analytic``
     - Local series, log branches, square roots and principal-unit powers
   * - ``extensions``
     - Certified unramified quadratic fields, common ideal-ball precision, exact quadratic residues and Frobenius
   * - ``ramified``
     - Eisenstein quadratic fields, weighted uniformizer bases, staggered coordinate precision and exact prime residues
   * - ``special``
     - Bernoulli values, theta series and explicit numerical KL branches
   * - ``symbolic``
     - Exact Q[x] and truncated Q[[x]], independent of optional CAS packages
   * - ``sympy_bridge``
     - Optional typed symbolic function expressions and safe numeric evaluation
   * - ``transforms``, ``measures``, ``taylor``
     - Exact measure data, bounded Amice transforms and certified local expansions
   * - ``harmonic``
     - Complex-place additive characters, real Haar integration and exact finite Fourier transforms
   * - ``wavelets``, ``operators``
     - Exact finite Kozyrev bases, whole-Qp derivatives with analytic tails, regional ball generators and numerical heat evolution
   * - ``algebra``
     - Scalar-generic quadratic, quaternion and diagonal Clifford algebras
   * - ``bridges``
     - Common exact sources, independent completion images and complex evaluation

The first four modules depend only on Python's standard library. Importing
``padic_tools`` does not import SymPy, mpmath or Sphinx. Exact rational
constants are handled before finite precision conversion wherever possible.
All extension points preserve the distinction between exact coefficients and
approximate p-adic balls.

The two quadratic field engines share ideal-ball arithmetic, with certified
basis weights and precision conversion supplied by their field classes.
Unramified bases are validated by irreducibility modulo p; ramified bases
by the Eisenstein criterion. Their valuations, residue degrees and norm/trace
precision are different. Read :doc:`extensions` and :doc:`ramified` before
adapting scalar-generic algebras or analytic functions to a larger field.

The harmonic layer also uses only the standard library for its exact
cyclotomic arithmetic. Its functions have rational/cyclotomic values in a
chosen complex place; Haar weights are real. Optional ``.complex(dps)``
evaluation uses an isolated mpmath context. Keep these integrals separate
from Qp-valued bounded measures and their Mellin/Amice transforms.

Adding an analytic function
---------------------------

Document the exact definition, scalar domain, branch choices, input and output
precision model, and a bound for the full omitted tail. State whether an
extension field is needed. Add tests against independent exact values or a
separately documented mathematical identity. A zero term or agreement between
two truncations of the same algorithm alone is not a tail proof.

Adding a bridge
---------------

State the common source algebra, each target algebra, choices of places or
generators, and which operations are preserved. Tests should compare
``image(a*b)`` with ``image(a)*image(b)`` and check inverses only when defined.
An analytic bridge also needs an argument that the relevant limits and
convergence domains are compatible. The quaternion/Clifford bridge is
algebraic and does not make that extra analytic claim.

Validation commands
-------------------

.. code-block:: console

   python -m pip install '.[dev]'
   python -m pytest
   python -m sphinx -W --keep-going -b html docs docs/_build/html
   python -m sphinx -W --keep-going -b doctest docs docs/_build/doctest
   python scripts/build_manual.py
   python -m build
   python -m twine check dist/*

The documented examples are executed by Sphinx doctest. CI also checks a
Python-version matrix, package builds and installation of the wheel with no
optional dependencies. The first release should additionally be compared
with Sage/PARI on matched conventions; the optional oracle test can run when
the ``gp`` executable is available.

The portable manual builder assembles the chapters into one file, embeds
styles and pre-renders formulas to native MathML using Node.js at build time.
It removes controls that depend on other files and validates every delivered
internal link. The resulting HTML can be read without JavaScript or network
access. Keep the complete directory together for the ordinary multi-page
edition; a link to only its index is not a portable manual.
