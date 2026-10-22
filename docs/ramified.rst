Ramified quadratic fields
==========================

``Q2i`` constructs Q2(i) in the uniformizer basis 1, pi, where pi=1+i.
``EisensteinQuadraticField`` provides the same arithmetic for any certified
quadratic Eisenstein polynomial over Qp, including wild dyadic cases.
Both use only the standard library.

.. doctest::

   >>> L = Q2i(precision=12)
   >>> pi = L.uniformizer()
   >>> i = pi-1
   >>> (i*i).is_congruent(-1)
   True
   >>> (pi*pi).is_congruent(2*i)
   True
   >>> (pi.valuation(), L(2).valuation(), pi.base_valuation())
   (1, 2, Fraction(1, 2))
   >>> pi.norm().contains(2), pi.trace().contains(2)
   (True, True)
   >>> (pi*pi.inverse()).is_congruent(1)
   True
   >>> L.residue_field().order
   2

The defining polynomial is

.. math::

   f(X)=X^2-2X+2,\qquad \pi=1+i,\qquad
   \pi^2=2i,\qquad e=2,\quad f_{\mathrm{res}}=1.

Its constant term has 2-adic valuation one and its linear coefficient is
divisible by 2, so it is Eisenstein. ``UnramifiedQuadraticField`` correctly
continues to reject X^2+1 at p=2; the new engine first selects a suitable
uniformizer basis. The two constructors certify different kinds of basis.

Defining data and valuation units
----------------------------------

For a general f(X)=X^2+B*X+C, the constructor requires exact rational
coefficients, monicity, v_p(B) at least one, and v_p(C)=1. These conditions
certify irreducibility and a totally ramified quadratic extension. The
generator pi is a uniformizer and the residue field is Fp. See `Sage's
Eisenstein extension reference
<https://doc.sagemath.org/html/en/reference/padics/sage/rings/padics/eisenstein_extension_generic.html>`_.

.. doctest::

   >>> R = EisensteinQuadraticField(Qp(3, 10), [-3, 0, 1], name="r")
   >>> (R.generator**2).is_congruent(3)
   True
   >>> (R.ramification_index, R.residue_degree)
   (2, 1)
   >>> EisensteinQuadraticField(Qp(2), [1, 0, 1])
   Traceback (most recent call last):
   ...
   padic_tools.core.DomainError: the polynomial must be Eisenstein: v_p(B)>=1 and v_p(C)=1

``valuation()`` uses v_pi(pi)=1, so v_pi(p)=2. The two power-basis terms
have different parity of valuation; they cannot cancel:

.. math::

   v_\pi(a+b\pi)=\min\{2v_p(a),\,2v_p(b)+1\},\qquad
   v_{\mathrm{base}}(x)=\tfrac12 v_\pi(x).

``base_valuation()`` returns the second normalization, with v_base(p)=1,
as an exact ``Fraction``. ``absolute_value_squared()`` returns p**(-v_pi),
the square of the absolute value extending the base p-adic norm. For odd
v_pi the absolute value itself is an irrational real square root; this API
returns its rational square without a floating-point rounding choice.
``norm()`` is instead the algebraic norm to Qp.

The Gaussian basis 1,i does not have these distinct valuation weights:
1+i has valuation one although both Gaussian coefficients are units. Supply
coordinates to ``L`` in the basis 1,pi, as in the final example below.

.. doctest::

   >>> (1+i).valuation(), (1-i).valuation()
   (1, 1)
   >>> pi.absolute_value_squared()
   Fraction(1, 2)
   >>> (1/pi).valuation()
   -1

Precision in powers of the uniformizer
---------------------------------------

An element is the entire ideal ball x0+pi^N O_L. Its ``precision`` and
``relative_precision`` count pi digits; ``absprec=N`` also uses pi units.
In the power basis this is exactly a rectangle of Qp coordinate balls:

.. math::

   a\pmod{p^{\lceil N/2\rceil}},\qquad
   b\pmod{p^{\lfloor N/2\rfloor}},\qquad
   N=\min\{2N_a,\,2N_b+1\}.

The last formula converts coordinate input precisions to the largest known
uniformizer ideal ball. It holds for negative N too. Odd precision retains
one additional base digit in the constant coordinate; rounding both
coordinates to the same precision would either lose information or invent it.

``Q2i(precision=P)`` takes a base Q2 precision P. Its field default is
2*P significant pi digits. More generally an exact input to the Eisenstein
field receives 2*base.precision significant pi digits. This scaling agrees
with `Magma's precision-of-extensions convention
<https://docs.magma-maths.org/LocalFields/PAdicRingsAndTheirExtensions/overview.html>`_.
``generator`` and ``uniformizer()`` therefore have absolute precision 2*P+1.
Existing Qp or extension balls cannot gain digits on conversion.

.. doctest::

   >>> (L.precision, pi.precision, pi.relative_precision)
   (24, 25, 24)
   >>> x = L(1, 1, absprec=5)
   >>> tuple(c.precision for c in x.coefficients)
   (3, 2)
   >>> limited = L(Qp(2)(1, absprec=4), Qp(2)(2, absprec=3))
   >>> limited.precision
   7
   >>> L(limited, absprec=8)
   Traceback (most recent call last):
   ...
   padic_tools.core.PrecisionError: cannot manufacture digits of an existing ball
   >>> z = L.O(-3)
   >>> tuple(c.precision for c in z.coefficients)
   (-1, -2)
   >>> z.contains(Fraction(1, 2), Fraction(1, 4))
   True

.. _ramified-precision-contract:

For input balls x and y with pi precisions N and M, and valuation lower
bounds v and w, arithmetic has the following guaranteed bounds. For a
zero ball its lower valuation bound is its absolute precision.

.. list-table:: Precision units are part of the contract
   :header-rows: 1
   :widths: 46 54

   * - Operation
     - Guaranteed absolute precision
   * - x+y, x-y
     - min(N,M) in pi units
   * - x*y
     - min(N+w,M+v) in pi units
   * - inverse(x), when zero is excluded
     - N-2*v in pi units
   * - conjugate(x)
     - N in pi units
   * - x times an exact nonzero q in Q
     - N+2*v_p(q) in pi units
   * - trace(x)
     - floor((N+d)/2) in base p units
   * - norm(x)
     - min(floor((N+v+d)/2),N) in base p units
   * - norm of a zero ball O(pi^N)
     - N in base p units

Norm and trace precision
-------------------------

Let d be the exponent of the different ideal. Since the two derivative
terms 2*pi and B have different parity of valuation, they do not cancel:

.. math::

   d=v_\pi(f'(\pi))=\min\{2v_p(2)+1,\,2v_p(B)\},\qquad
   \operatorname{Tr}(\pi^N\mathcal O_L)
      =p^{\lfloor(N+d)/2\rfloor}\mathbb Z_p.

The trace-ideal formula can be checked directly in our basis. The trace of
a+b*pi is 2*a-B*b. Its two coefficient errors have base valuation bounds
v_p(2)+ceil(N/2) and v_p(B)+floor(N/2). Their minimum is floor((N+d)/2),
and independent coefficient choices attain that trace-ideal bound.

For a perturbation delta in pi^N O_L, the norm difference is
trace(sigma(x)*delta)+norm(delta). Multiplication by sigma(x) carries the
first error into pi^(N+v) O_L. The trace-ideal bound gives its base
precision floor((N+v+d)/2). The second error has base valuation at least N,
because the residue degree is one and v_p(norm(delta))=v_pi(delta).
Taking the minimum proves the norm bound. It is a guaranteed enclosure;
dyadic cancellations can sometimes give more digits than this bound.
For a zero ball we obtain O(p^N).

Q2(i) has d=2; Q2(sqrt(2)) in its Eisenstein basis has d=3. The code
handles both wild cases, and the tame odd-prime case d=1.

.. doctest::

   >>> (L.different_exponent, EisensteinQuadraticField(Qp(2), [-2, 0, 1]).different_exponent)
   (2, 3)
   >>> x = L(0, 2, absprec=11)
   >>> (x.valuation(), x.inverse().precision, (x*2).precision)
   (3, 5, 13)
   >>> (x.trace().precision, x.norm().precision)
   (6, 8)
   >>> (L.O(-3).trace().precision, L.O(-3).norm().precision)
   (-1, -3)
   >>> L.O(2).inverse()
   Traceback (most recent call last):
   ...
   padic_tools.core.PrecisionError: cannot invert an extension ball containing zero

Residues, conjugation and the Gaussian source
----------------------------------------------

An integral element reduces to a modulo p, since pi reduces to zero. One
known pi digit suffices: the constant coordinate then has at least one
known base digit. ``PrimeResidueField`` and ``PrimeResidueElement`` implement
exact arithmetic, inversion and Frobenius on Fp. Integer coefficient lifting
is a section of reduction and does not preserve all ring operations.

.. doctest::

   >>> F2 = L.residue_field()
   >>> L(3, 1, absprec=1).residue() == F2(1)
   True
   >>> pi.residue() == F2(0)
   True
   >>> F2(1)**2 == F2(1).frobenius()
   True
   >>> L.lift_residue(F2(1)).residue() == F2(1)
   True
   >>> pi.conjugate().is_congruent(2-pi)
   True

The local-field conjugation sends pi to -B-pi. Its action on the prime
residue field is the identity, as is residue Frobenius. There is no unique
local Frobenius lift singled out by that residue action in a ramified
extension, so the local element has ``conjugate()`` and no ``frobenius()``
method. The unramified field retains its canonical Frobenius API.

For an exact Gaussian source a+b*i, choose i=pi-1 and change the coordinates
to (a-b,b). This defines its local image without using a complex
approximation or a finite residue as exact source data:

.. doctest::

   >>> source = Gaussian(Fraction(1, 2), 3)
   >>> image = L(source.a-source.b, source.b)
   >>> image.norm().contains(source.norm())
   True
   >>> square = source*source
   >>> (image*image).contains(square.a-square.b, square.b)
   True

Dedicated Gaussian/quaternion bridge adapters are the next step. The
existing coefficient-algebra bridge methods keep their documented return
types. Arbitrary ramified polynomials need a certified change to an
integral uniformizer basis; higher extension degrees, extension analytic
functions and matrix solvers also remain subsequent work.

Validation covers odd and negative pi precision, whole-ideal perturbations,
exact small quotient sizes, independent SymPy polynomial inverses and
multiplication-matrix norm/trace, regression tests of the unramified engine,
and PARI/GP polymod checks at even and odd uniformizer valuations. The PARI
valuation oracle uses the norm to Qp and residue degree one, rather than
applying a coefficient-minimum valuation to a ramified polynomial basis.
