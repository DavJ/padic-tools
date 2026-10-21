Numerical tutorial
==================

Construction and precision
--------------------------

.. doctest::

   >>> K = Qp(5, 12)
   >>> a = K(25)
   >>> (a.valuation(), a.precision, a.relative_precision)
   (2, 14, 12)
   >>> b = K("1/25")
   >>> (b.valuation(), b.precision)
   (-2, 10)
   >>> a.norm()
   Fraction(1, 25)
   >>> K(Fraction(1, 3), absprec=6).digits()
   (0, (2, 3, 1, 3, 1, 3))

Use integers, ``Fraction``, ``Decimal`` or text. Implicit binary floats and
Python ``complex`` are rejected. If a binary float really is intended,
explicitly pass ``Fraction.from_float(value)``.

``absprec`` sets an absolute precision for a fresh exact input. Applying it
to an existing ball can only discard information. A zero ball has a valuation
lower bound, not a known infinite valuation.

Comparison
----------

.. doctest::

   >>> a, b = K(1, absprec=2), K(26, absprec=5)
   >>> a == b
   False
   >>> a.is_congruent(b)
   True
   >>> a.contains(26)
   True
   >>> K.O(7).is_zero()
   True

Equality is equality of balls, with transitive semantics and consistent
hashes. ``is_congruent`` compares known digits; ``contains`` tests an exact
rational. A truth-value test such as ``if x`` raises because it would conceal
this distinction. The rational returned by ``representative()`` is only
one possible centre, not a recovered exact number.

Analytic functions
------------------

.. doctest::

   >>> x = K(5)
   >>> log(exp(x)).is_congruent(x)
   True
   >>> (sin(x)**2 + cos(x)**2).is_congruent(1)
   True
   >>> power(1+x, 3).is_congruent((1+x)**3)
   True
   >>> root = sqrt(K(9))
   >>> (root*root).is_congruent(9)
   True
   >>> log(K(5), p_branch=7).is_congruent(7)
   True

``sqrt`` chooses one documented residue branch. It need not return the
positive real root: positivity is an archimedean concept. Requesting more
digits than are justified raises ``PrecisionError``.

Hensel lifting
--------------

.. doctest::

   >>> f = Polynomial([-6, 0, 1])  # x^2-6, with root 1 modulo 5
   >>> r = hensel_lift(f, 1, K)
   >>> f(r).is_zero()
   True

The native routine implements the simple-root theorem. It rejects singular
residue roots and coefficients outside Zp rather than applying the wrong
lifting theorem. Polynomial coefficients remain exact rational numbers.
