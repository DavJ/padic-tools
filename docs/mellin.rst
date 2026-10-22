Mellin transforms of rational measures
=======================================

``RationalAmiceMeasure`` supports the same character convention as
``AtomicMeasure`` in :doc:`transforms`, now for a complete rational
Amice transform. Restriction to units is exact, and numerical integration
includes a proved bound for the entire omitted Mahler tail.

Exact restriction to units
---------------------------

Write z=1+T and F(z)=A_mu(z-1). The root-of-unity filter gives

.. math::

   A_{\mu|\mathbb Z_p^\times}(T)
   = A_\mu(T)-\frac1p\sum_{\zeta^p=1}A_\mu(\zeta(1+T)-1).

The average selects the Laurent powers z^n with n divisible by p,
including negative powers. This identity is the unit-restriction formula
in Rodrigues Jacinto and Williams, section 2.4; see :doc:`references`.

The implementation computes that rational section over Q. After removing
a possible finite pole at z=0 and a polynomial part, the remaining
coefficients satisfy a finite linear recurrence. If C is its companion
matrix, stepping by p replaces C by C^p. The characteristic polynomial of
C^p recovers the exact decimated rational function. Subtracting this section
and replacing z by 1+T produces another checked rational Amice measure.
No numerical roots of unity or truncated power-series inference are used.

.. doctest::

   >>> K = Qp(5, 18)
   >>> mu = RationalAmiceMeasure(K, 1, [2, 1])
   >>> units = mu.restrict_units()
   >>> units.moment(0), units.moment(1), units.moment(3)
   (Fraction(0, 1), Fraction(1, 1), Fraction(-31, 2))
   >>> units.restrict_units() is units
   True
   >>> atoms = AtomicMeasure(K, [(-2, 1), (0, 7), (5, 3), (6, 2)])
   >>> RationalAmiceMeasure.from_atomic(atoms).restrict_units() == RationalAmiceMeasure.from_atomic(atoms.restrict_units())
   True

``restrict_units(max_degree=4096)`` limits the estimated degree of section
polynomials; the native recurrence order is at most 64. Exceeding these
resource limits raises ``PrecisionError``. A restriction already proved to
have unit support returns immediately. Addition and scaling preserve that
proof; additive convolution requires checking again because sums of units
can be nonunits.

Character and domain
---------------------

.. math::

   M_{\mu,r}(s)=\int_{\mathbb Z_p^\times}
      \omega(x)^r\langle x\rangle^s\,d\mu(x),\qquad s\in\mathbb Z_p.

At odd primes r is reduced modulo p-1. At p=2 the character omega is the
sign modulo 4, and r is reduced modulo 2. The exponent has no hidden shift.
Nonunit support raises ``DomainError``; call ``restrict_units`` explicitly
when that change of measure is intended. For a nonnegative integer n in
branch n, the transform interpolates the restricted exact n-th moment.

.. doctest::

   >>> units.mellin(3, branch=3).contains(units.moment(3))
   True
   >>> units.mellin(Fraction(1, 2), branch=3).precision
   18
   >>> dyadic = RationalAmiceMeasure(Qp(2, 14), 1, [3, 1]).restrict_units()
   >>> dyadic.mellin(3, branch=1).contains(dyadic.moment(3))
   True

Mahler tail certificate
-------------------------

Extend the kernel by zero on nonunits and expand it in binomial polynomials:

.. math::

   f_{s,r}(x)=\sum_{n\geq0}c_n{x\choose n},\qquad
   v_p(c_n)\geq g_p(n),\qquad
   g_p(n)=\begin{cases}v_p(n!),&p\text{ odd},\\
                       \lfloor n/2\rfloor-1,&p=2.
          \end{cases}

Roblot's Proposition 5.10 and Corollary 5.13 prove these bounds for the
principal branch. The same proof covers our torsion branches: at odd p,
choose arbitrarily large nonnegative integers t tending to s p-adically
with t congruent to r modulo p-1. Their polynomials x^t tend uniformly to
the zero-extended kernel. Their n-th differences are multiples of n!,
so the bound passes to the limit. At p=2, matching parity gives polynomial
approximants; the other parity multiplies them by the sign character.
Its k-th Mahler coefficient has valuation at least floor(k/2).
Roblot's product-coefficient lemma gives floor(n/2)-1. These arguments
apply to all s in Zp and all implemented branches.

Let b be ``valuation_bound`` for the measure. To return absolute precision
P, set M=P-b. If M is positive, choose N with g_p(N)>=M. Monotonicity
bounds every omitted coefficient, hence

.. math::

   M_{\mu,r}(s)-\sum_{n=0}^{N-1}c_na_n\ \in O(p^{M+b})=O(p^P).

For odd p the evaluator finds the least N with v_p(N!)>=M; for p=2 it
uses N=2(M+1). Kernel samples modulo p^M give Mahler coefficients by
forward differences, then are paired with the exact Amice coefficients
a_n. Modular integer powers evaluate the character without repeatedly
computing logarithms and exponentials. This is a finite algorithm with a
bound for its infinite tail, rather than a numerical stabilization test.

Precision and work limits
---------------------------

``precision=P`` means absolute output order, including for measures with
nonintegral weights. It defaults to the field's precision. If an exponent
is already a ball s+O(p^S), its uncertainty permits output order at most
S+e+b, with e=1 at odd primes and e=2 at p=2. This follows from
``v_p(log(<x>))>=e``. The default is capped accordingly; explicitly
requesting more raises ``PrecisionError``. Exact rational exponents can
instead be evaluated with additional working digits.

.. doctest::

   >>> low = K(Fraction(1, 2), absprec=3)
   >>> units.mellin(low, branch=3).precision
   4
   >>> units.mellin(Fraction(1, 2), branch=3, precision=22).precision
   22
   >>> units.mellin(low, branch=3, precision=8)
   Traceback (most recent call last):
       ...
   padic_tools.core.PrecisionError: Mellin precision exceeds the available exponent digits

``max_terms=10000`` caps the number N of retained Mahler coefficients.
Insufficient terms raise ``PrecisionError``; the evaluator never returns
an uncertified cutoff. The native forward-difference stage uses quadratic
work in N, so large primes or requested precisions can be expensive.
Increasing the cap is an explicit resource choice, not a convergence proof.

A check against Kubota--Leopoldt zeta
-------------------------------------

For odd p, the measure A(T)=1/(2+T) in :doc:`measures` satisfies
``m_k=(-1)^k*(1-2**(k+1))*zeta_special(-k)``. Restricting it to units gives
the exact Euler factor ``1-p**k``. Consequently, on branch r=p-2,

.. math::

   M_{\mu|\mathbb Z_p^\times,p-2}(s)
   =\bigl(\langle2\rangle^{s+1}-1\bigr)
        \zeta_{p,0}(-s),\qquad s\ne-1.

To check the shifts, take positive n congruent to p-2 modulo p-1.
Then n is odd, n+1 is divisible by p-1, and our principal KL branch at
-n is ``(1-p**n)*zeta_special(-n)``. The moment identity gives the displayed
formula; continuity extends it from this dense integer set. The expression
at s=-1 has a removable singularity, while direct Mellin integration
remains defined there. It does not call ``zeta_p`` at its pole.

.. doctest::

   >>> s = Fraction(1, 2)
   >>> principal2 = K(2)/teichmuller(K(2))
   >>> expected = (power(principal2, s+1)-1)*zeta_p(-s, K, branch=0)
   >>> units.mellin(s, branch=3).is_congruent(expected)
   True
   >>> units.mellin(-1, branch=3).precision
   18

The tests additionally compare this noninteger value to PARI/GP's
independent KL implementation at p=3, 5 and 7. Exact unit restriction is
also checked by inversion in Q(z)[w]/(w^p-1), independently of companion
recurrences. These are checks within the stated p-adic conventions;
they do not supply a correspondence between complex and p-adic zeros.
