Special functions: definitions and algorithms
=============================================

Elementary functions
--------------------

``exp``, ``sin`` and ``cos`` use their rational Taylor coefficients on
``v_p(x)>1/(p-1)``. ``log1p(x)`` requires ``v_p(x)>=1``. General ``log(x)``
removes the p-power and root-of-unity factors and applies the local series,
with the explicit branch ``log(p)=p_branch``. At p=2 the unit sign is chosen
modulo 4. ``power`` accepts only principal units in the exponential domain
and exponents in Zp. The tail arguments are in :doc:`precision`; compare
the `Sage analytic-function documentation <https://doc.sagemath.org/html/en/reference/padics/sage/rings/padics/padic_generic_element.html>`_.

Multiplicative theta
--------------------

The native two-argument definition is

.. math::

   \Theta(q,z)=\sum_{n\in\mathbb Z}q^{n^2}z^n
   =1+\sum_{n\ge1}q^{n^2}(z^n+z^{-n}).

For ``v_p(q)>0`` this converges for every nonzero z in Qp. In the complex
domain it agrees with Jacobi theta_3(w,q) after the explicit substitution
``z=exp(2*i*w)``; see `DLMF 20.2.3 <https://dlmf.nist.gov/20.2.E3>`_.
That complex substitution is not silently imposed p-adically.
``theta4(q)`` is ``Theta(q,-1)``.

Let a=v_p(q)>0 and b=v_p(z). Every omitted pair at index n has valuation
at least ``a*n*n-|b|*n``. The evaluator stops only when this bound has
reached the target and is nondecreasing for all later n. This includes z
with nonzero valuation. Default precision also accounts for amplification
of input errors by these nonintegral z powers.

Only theta3 and theta4 in this stated convention are implemented natively.
Theta1/theta2 require explicit fractional nome powers or an extension and
branch convention. Multidimensional/rigid-analytic theta functions are
separate future features, not aliases for this function.

Three distinct zeta APIs
------------------------

.. list-table::
   :header-rows: 1
   :widths: 25 35 40

   * - API
     - Domain / codomain
     - Meaning
   * - ``zeta_special(s)``
     - Integers s<=0 -> exact Q
     - Classical Riemann zeta special value
   * - ``riemann_zeta(s, dps=...)``
     - Exact Q or Q[i] input -> numerical C
     - Classical complex zeta via optional mpmath
   * - ``zeta_p(s, K, branch=r)``
     - Zp away from s=1 -> finite precision Qp
     - The explicitly defined even KL branch below

With B_1=-1/2, ``zeta_special(-m)=(-1)^m*B_(m+1)/(m+1)`` for m>=0.
This agrees with `DLMF 25.6 <https://dlmf.nist.gov/25.6>`_, including the
separate value zeta(0)=-1/2. ``zeta_euler_value`` multiplies by
``1-p**(-s)`` but does not by itself choose an analytic branch.

KL branch convention
--------------------

For odd p write ``a=omega(a)*<a>`` using the Teichmüller unit. Put F=p.
For p=2 put F=4, use a in {1,3}, and use omega(a)=+1 or -1 according to
the residue modulo 4. For even r modulo p-1 (only r=0 at p=2), define

.. math::

   Z_{p,r}(s)=\frac{1}{F(s-1)}
   \sum_{\substack{1\le a\lt F\\p\nmid a}}
   \omega(a)^r\langle a\rangle^{1-s}
   \sum_{k\ge0}\binom{1-s}{k}\left(\frac Fa\right)^k B_k.

This follows by expanding the Volkenborn integral on each residue disc;
``integral(t^k dt)=B_k``. The integral identities and Teichmüller
decomposition are recalled in `Lai, Lupu and Sprang, sections 2.1–2.3
<https://link.springer.com/article/10.1007/s40687-025-00559-x>`_. Our displayed
formula, residue-disc normalization and implementation convention are
spelled out here to make comparisons reproducible.

At s=1-n with n>=2 and n congruent to r modulo p-1 (even n for p=2),
``omega(a)^r*<a>^n=a^n``. The finite Bernoulli-polynomial multiplication
formula then gives

.. math::

   Z_{p,r}(1-n)=(1-p^{n-1})\zeta(1-n).

This is the special-value bridge; see also Rodrigues Jacinto and Williams,
`Theorem 3.12 and section 4.3 <https://warwick.ac.uk/fac/sci/maths/people/staff/cwilliams/lecturenotes/lecture_notes_part_i.pdf>`_.
The package's ``branch=0`` fixes r while s varies. Some literature instead
uses a character depending on an integer argument s. Those conventions can
give different values at positive integers; the branch formula above is the
contract. Match conventions before comparing to PARI, Sage or a paper.

.. doctest::

   >>> K = Qp(5, 16)
   >>> zeta_p(-1, K, branch=2).is_congruent(Fraction(1, 3))
   True
   >>> zeta_p(-3, K, branch=0).is_congruent(Fraction(-31, 30))
   True
   >>> zeta_p(-1, K, branch=0).is_congruent(zeta_euler_value(-1, 5))
   False
   >>> zeta_p(Fraction(1, 6), K).precision
   16

Tail certification
------------------

For s in Zp, ``binomial(1-s,k)`` is p-integral, and von Staudt–Clausen gives
``v_p(B_k)>=-1``. Let e=v_p(F) and d=v_p(s-1). The entire tail starting at
k=K therefore has valuation at least

.. math::

   Ke-1-e-d.

The implementation chooses K so this is at least the requested output
precision. The finite sum uses ordinary ball arithmetic, tracking all
divisions and the pole denominator separately. Its result is then truncated
to the proven tail precision. Exact rational input permits guard digits;
an existing finite p-adic input does not. Near s=1 the latter can lose
substantial precision.

Pole-containing balls are rejected. At s=1 branch 0 has a pole;
nonprincipal even branches have a removable limit whose evaluation is not
implemented. General Dirichlet characters and p-adic arguments in extension
fields are outside 0.1. Complexity scales with the number of residue discs
and target precision, so this native routine is intended for modest primes.

The ordinary Dirichlet series ``sum(n**(-s))`` is not used: its summands
generally do not tend to zero in the p-adic norm. No complex functional
equation or correspondence between complex and p-adic zeros is assumed.
