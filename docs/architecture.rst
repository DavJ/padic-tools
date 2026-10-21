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
   * - ``special``
     - Bernoulli values, theta series and explicit numerical KL branches
   * - ``symbolic``
     - Exact Q[x] and truncated Q[[x]], independent of optional CAS packages
   * - ``sympy_bridge``
     - Optional typed symbolic function expressions and safe numeric evaluation
   * - ``algebra``
     - Scalar-generic quadratic, quaternion and diagonal Clifford algebras
   * - ``bridges``
     - Common exact sources, independent completion images and complex evaluation

The first four modules depend only on Python's standard library. Importing
``padic_tools`` does not import SymPy, mpmath or Sphinx. Exact rational
constants are handled before finite precision conversion wherever possible.
All extension points preserve the distinction between exact coefficients and
approximate p-adic balls.

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
   python -m build
   python -m twine check dist/*

The documented examples are executed by Sphinx doctest. CI also checks a
Python-version matrix, package builds and installation of the wheel with no
optional dependencies. The first release should additionally be compared
with Sage/PARI on matched conventions; the optional oracle test can run when
the ``gp`` executable is available.
