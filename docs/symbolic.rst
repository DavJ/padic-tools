Symbolic polynomials, series and functions
==========================================

Native exact layer
------------------

``Polynomial`` represents Q[x], with ascending rational coefficients.
It provides exact arithmetic, Euclidean division, monic gcd, differentiation,
integration and evaluation in a scalar algebra. ``PowerSeries`` represents
Q[[x]] modulo x^N, with its own formal order separate from numerical p-adic
precision.

``Polynomial.shift(center)`` translates a polynomial exactly. ``taylor``
adds rational-centered jets and certified local numerical expansions;
see :doc:`taylor`. Exact moment series of measures are described in
:doc:`transforms`.

.. doctest::

   >>> x = Polynomial([0, 1])
   >>> f = (x-1)**3*(x+2)
   >>> f.gcd(f.derivative()) == (x-1)**2
   True
   >>> t = PowerSeries([0, 1], order=12)
   >>> t.exp().log() == t
   True
   >>> (1+t).log() == log1p_series(12)
   True
   >>> theta3_series(10).coefficients
   (Fraction(1, 1), Fraction(2, 1), Fraction(0, 1), Fraction(0, 1), Fraction(2, 1), Fraction(0, 1), Fraction(0, 1), Fraction(0, 1), Fraction(0, 1), Fraction(2, 1))

Formal exp requires constant coefficient zero, and formal log requires
constant one. Composition requires inner constant zero. A derivative loses
one known coefficient; an integral gains one. ``polynomial()`` explicitly
chooses a finite truncation and does not promise that it approximates an
infinite function at a numerical argument. Use the analytic APIs when a
tail enclosure is needed.

Optional SymPy layer
--------------------

.. doctest::

   >>> import sympy as sp
   >>> from padic_tools.sympy_bridge import PadicExp, PadicLog, PadicTheta3, PadicZeta, evaluate
   >>> s = sp.Symbol("s")
   >>> sp.diff(PadicExp(s), s) == PadicExp(s)
   True
   >>> sp.diff(PadicLog(s), s)
   1/s
   >>> evaluate(PadicExp(s), Qp(5), {s: 5}).is_congruent(exp(Qp(5)(5)))
   True
   >>> evaluate(PadicZeta(-3, 5, 0), Qp(5)).is_congruent(Fraction(-31, 30))
   True

The typed function nodes keep domain information in their names and prime /
branch data in zeta's arguments. SymPy provides exact rearrangement and
differentiation; the adapter recursively evaluates supported expressions
with the numerical library's precision rules and domain checks.
It rejects ordinary SymPy ``zeta``, unbound symbols, implicit complex
constants, floats and noninteger powers. A general symbolic zero test,
algebraic closure or theorem prover is not claimed.

Quaternion and Clifford formulas can also be symbolic by choosing
``QuaternionAlgebra(base=sp.sympify)`` or
``CliffordAlgebra((1,1,1), base=sp.sympify)``. The scalar expressions should
commute; base-dependent simplification is the caller's responsibility.
Numerical evaluation is an explicit base change, never an implicit coercion
of a complex approximation into Qp. See the `SymPy custom-function guide
<https://docs.sympy.org/latest/guides/custom-functions.html>`_ for the adapter
mechanism.
