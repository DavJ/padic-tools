Bounded p-adic measures
=======================

Rational Amice data
--------------------

``RationalAmiceMeasure`` specifies a bounded Qp-valued measure on Zp by
its complete rational Amice transform, rather than a finite atom list:

.. math::

   A_\mu(T)=\frac{P(T)}{Q(T)}=\sum_{n\geq0}a_nT^n,
   \qquad a_n=\int_{\mathbb Z_p}{x\choose n}\,d\mu(x).

Numerator and denominator are exact rational polynomials. Common factors
are cancelled before domain checks, then Q is normalized to Q(0)=1.
Every other normalized denominator coefficient must be p-integral.
The constructor verifies this condition; it does not infer boundedness
from a finite numerical sample.

If b is the minimum numerator coefficient valuation, the recurrence
``a_n=P_n-sum(Q_j*a_(n-j), j>=1)`` proves ``v_p(a_n)>=b`` by induction.
The measure/Amice correspondence then gives a bounded measure, after
scaling to the integral-coefficient case. Zero measures use the safe
bound b=0. The bound is sufficient and need not be optimal.

On v_p(T)>=1, Q(T) is a unit: its constant coefficient is one and every
other term has positive valuation. Numerical evaluation uses the complete
rational expression, so it has no artificial series cutoff.

.. doctest::

   >>> K = Qp(5, 18)
   >>> mu = RationalAmiceMeasure(K, [1], [2, 1])
   >>> mu.amice_series(5).coefficients
   (Fraction(1, 2), Fraction(-1, 4), Fraction(1, 8), Fraction(-1, 16), Fraction(1, 32))
   >>> mu.valuation_bound
   0
   >>> mu.amice(K(5)).is_congruent(Fraction(1, 7))
   True
   >>> mu.laplace(K(5)).is_congruent(1/(1+exp(K(-5))))
   True

The example ``1/(2+T)`` passes the boundedness criterion for odd p.
It fails at p=2 because its coefficients have unbounded negative
2-adic valuation. For instance ``1/(3+T)`` is available at p=2.
Scalar denominator powers of p are allowed when normalization only
rescales the numerator; ``1/p`` is a bounded constant transform.

Moments, integration and convolution
--------------------------------------

Ordinary moments are recovered exactly using Stirling numbers of the
second kind, denoted S(n,k):

.. math::

   \int x^n\,d\mu(x)=\sum_{k=0}^n S(n,k)\,k!\,a_k.

Polynomial integration follows by linearity. Additive convolution
multiplies Amice transforms; addition and exact rational scaling are
also implemented. Operations check that the prime agrees and preserve
the bounded rational class.

.. doctest::

   >>> mu.moment(0), mu.moment(1), mu.moment(2), mu.moment(3)
   (Fraction(1, 2), Fraction(-1, 4), Fraction(0, 1), Fraction(1, 8))
   >>> mu.integrate_polynomial(Polynomial([1, 2, 3]))
   Fraction(0, 1)
   >>> mu.convolve(mu).laplace(K(5)).is_congruent(mu.laplace(K(5))**2)
   True
   >>> mu.scale(2).moment(0)
   Fraction(1, 1)
   >>> atom = AtomicMeasure(K, [(-2, 1), (3, 2)])
   >>> converted = RationalAmiceMeasure.from_atomic(atom)
   >>> converted.amice_series(12) == atom.amice_series(12)
   True

Conversion from atomic data is supported for integer points, including
negative integers through denominator powers of 1+T. Noninteger rational
points generally require algebraic transforms and are explicitly rejected.

Mahler coefficients with an explicit tail contract
---------------------------------------------------

A continuous function on Zp has a Mahler expansion in binomial polynomials.
For known exact coefficients c_0,...,c_(N-1), suppose the caller has
established a uniform lower bound B on the valuations of all omitted c_n.
Then

.. math::

   \int f\,d\mu-\sum_{n=0}^{N-1}c_na_n\ \in O(p^{B+b}).

``integrate_mahler`` returns this enclosure. The caller must establish
continuity and the infinite tail bound; the library cannot prove them
from a finite coefficient prefix. An omitted tail is never silently
assumed to be zero. ``integrate_polynomial`` handles exact finite input
without that additional contract.

.. doctest::

   >>> # (1+5)^x has exact Mahler coefficients 5^n, with tail valuation >=N.
   >>> value = mu.integrate_mahler([5**n for n in range(12)], tail_valuation=12)
   >>> value.precision
   12
   >>> value.is_congruent(mu.amice(5), precision=12)
   True

An exact generating-function bridge
------------------------------------

For odd p and A(T)=1/(2+T), the common rational generating function gives
the following exact rational moment identity for every integer k>=0:

.. math::

   m_k=(-1)^k(1-2^{k+1})\zeta(-k).

Indeed the exponential generating series of these moments is
``A(exp(t)-1)=1/(1+exp(t))``; the Bernoulli generating formula gives
the displayed coefficients. The local Laplace convention instead uses
``A(exp(-t)-1)``. Thus its Taylor coefficients are ``(-1)^k*m_k/k!``.

.. doctest::

   >>> all(mu.moment(k) == (-1)**k*(1-2**(k+1))*zeta_special(-k) for k in range(15))
   True

This compares exact moment data of the rational generating function.
It is not a theorem identifying zeros of complex and p-adic zeta functions.
See the measure/Amice references in :doc:`references` for the underlying
correspondence, and the Bernoulli/zeta source there for the special values.

Scope
------

This class implements bounded rational Amice measures, not every bounded
measure on Zp. General Mellin evaluation and restriction to units for this
class, arbitrary coefficient oracles, and inverse transforms remain future
work. ``AtomicMeasure.mellin`` continues to support finite unit measures.
The local Laplace domain here is the sufficient domain for all of Zp:
``v_p(t)>1/(p-1)``. Special support can justify wider domains in other APIs.
