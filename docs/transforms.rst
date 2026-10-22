Mellin, Laplace and Amice transforms
====================================

An explicit measure
-------------------

``AtomicMeasure`` represents a finite, bounded Qp-valued measure on Zp:

.. math::

   \mu=\sum_{j=1}^m w_j\delta_{a_j},\qquad
   \int f(x)\,d\mu(x)=\sum_{j=1}^m w_j f(a_j).

Points and weights are exact rational source data. Points must be
p-integral; weights may be negative or nonintegral. Repeated points are
combined and zero weights removed. This is a concrete measure, rather than
an unspecified integration rule or an approximation to an arbitrary measure.
It is not automatically a probability measure.

.. doctest::

   >>> K = Qp(5, 18)
   >>> mu = AtomicMeasure(K, [(1, 1), (2, 2)])
   >>> mu.moment(0), mu.moment(1), mu.moment(2)
   (Fraction(3, 1), Fraction(5, 1), Fraction(9, 1))

Local p-adic Laplace transform
------------------------------

Our sign convention is

.. math::

   L_\mu(t)=\int_{\mathbb Z_p}\exp_p(-tx)\,d\mu(x),\qquad
   L_\mu(t)=\sum_{n\geq0}\frac{(-1)^n m_n}{n!}t^n,
   \quad m_n=\int x^n\,d\mu(x).

Every nonzero atom must satisfy ``v_p(t*a_j)>1/(p-1)``. In particular,
``t`` in pZp suffices at odd p, and ``t`` in 4Z2 suffices at p=2.
Support in a smaller ball can allow more parameters. Each kernel evaluation
uses the certified exponential implementation; exact weight multiplication
can reduce the output precision. An atom at zero contributes its exact weight.

.. doctest::

   >>> t = K(5)
   >>> mu.laplace(t).is_congruent(exp(-t) + 2*exp(-2*t))
   True
   >>> mu.laplace_series(4).coefficients
   (Fraction(3, 1), Fraction(-5, 1), Fraction(9, 2), Fraction(-17, 6))
   >>> AtomicMeasure(Qp(2), [(4, 1)]).laplace(1).is_congruent(exp(Qp(2)(-4)))
   True

``laplace_series`` returns exact formal coefficients modulo t^order.
Numerical truncation of that polynomial alone does not carry an analytic
remainder. Use ``laplace`` for the enclosed value. There is no assumed
ordered half-line in Qp, nor an implicit complex Bromwich inversion.

Branched p-adic Mellin transform
--------------------------------

For a measure supported on units, use the explicit convention

.. math::

   M_{\mu,r}(s)=\int_{\mathbb Z_p^\times}
      \omega(x)^r\langle x\rangle^s\,d\mu(x),\qquad
   \langle x\rangle=x/\omega(x),\quad s\in\mathbb Z_p.

At odd p, omega is the Teichmuller character and r is taken modulo p-1.
At p=2 it is the sign determined modulo 4, with r modulo 2; the principal
part lies in 1+4Z2. Principal-unit powers are evaluated as exp(s*log(<x>)).
The exponent s is used as written, with no hidden s-1 or 1-s shift.

For a nonnegative integer n congruent to r modulo the branch period,
``mellin(n, branch=r)`` equals the exact n-th moment. Noninteger p-adic
exponents are also supported.

.. doctest::

   >>> mu.mellin(6, branch=2).is_congruent(mu.moment(6))
   True
   >>> mu.mellin(Fraction(1, 2), branch=2).precision >= 18
   True
   >>> units2 = AtomicMeasure(Qp(2, 18), [(-1, 1), (5, 2)])
   >>> units2.mellin(3, branch=1).is_congruent(units2.moment(3))
   True
   >>> mixed = AtomicMeasure(K, [(0, 7), (5, 2), (1, 4)])
   >>> mixed.restrict_units().mellin(0).is_congruent(4)
   True

``mellin`` rejects nonunit support. Restriction is explicit because it
changes the measure and its moments. The branch convention here describes
the character of a general measure; it is not the argument convention of
``zeta_p``. No correspondence of complex and p-adic zeros is asserted.

Amice transform and its Laplace relation
----------------------------------------

.. math::

   A_\mu(T)=\int(1+T)^x\,d\mu(x)
     =\sum_{n\geq0}\left(\int {x\choose n}\,d\mu(x)\right)T^n,
   \qquad L_\mu(t)=A_\mu(\exp_p(-t)-1).

The last identity holds formally and on the common numerical domain.
The binomial series converges on ``v_p(T)>=1``, including T in 2Z2.
This domain is larger than the domain of the dyadic exponential.

If every weight has valuation at least b, then every Amice coefficient
has valuation at least b: binomial(x,n) is p-integral for x in Zp.
The full tail from order N therefore belongs to
``O(p**(b+N*v_p(T)))``. The numerical evaluator adds this enclosure to
the finite polynomial evaluation and also propagates input uncertainty.

.. doctest::

   >>> u = PowerSeries([0, 1], order=12)
   >>> mu.amice_series(12).compose((-u).exp()-1) == mu.laplace_series(12)
   True
   >>> mu.amice(exp(-t)-1).is_congruent(mu.laplace(t))
   True
   >>> third = AtomicMeasure(Qp(2, 18), [(Fraction(1, 3), 1)])
   >>> (third.amice(2)**3).is_congruent(3)
   True

Scope and sources
-----------------

``AtomicMeasure`` covers finite atomic measures with rational source data.
``RationalAmiceMeasure`` additionally supports measures specified by a
complete rational Amice transform, exact moments, local Laplace evaluation,
exact unit restriction and Mellin evaluation with proved Mahler tails; see
:doc:`measures` and :doc:`mellin`. General bounded measures and transform inverses remain
future extensions. Qp-valued bounded measures and the real Haar measure
used for complex-valued harmonic analysis are separate integration settings;
the latter is implemented for exact finite test functions in :doc:`harmonic`.

For the measure/series correspondence, Dirac measures and Mellin characters,
see Rodrigues Jacinto and Williams, sections 2.2--2.4 and 4.3, linked in
:doc:`references`. The Amice coefficient identity is also documented in
Mathlib's ``NumberTheory.Padics.Measure.AmiceTransform``. The bounds above
specialize these definitions to explicit finite measures.
