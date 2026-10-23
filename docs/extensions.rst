Unramified quadratic fields
===========================

``UnramifiedQuadraticField`` constructs a genuine quadratic extension of Qp
with a certified integral power basis, an ideal-ball precision model and
the exact residue field of order p squared. It uses only the standard library.
Start with Q3(i):

.. doctest::

   >>> L = UnramifiedQuadraticField(Qp(3, 12), [1, 0, 1], name="i")
   >>> i = L.generator
   >>> (i*i).is_congruent(-1)
   True
   >>> (L.degree, L.ramification_index, L.residue_degree)
   (2, 1, 2)
   >>> L.uniformizer().valuation()
   1
   >>> x = L(1, 2)
   >>> (x*x.inverse()).is_congruent(1)
   True
   >>> x.norm().contains(5), x.trace().contains(2)
   (True, True)

The same polynomial is split over Q5 and ramified over Q2. Neither is an
unramified quadratic presentation, so this constructor rejects both:

.. doctest::

   >>> UnramifiedQuadraticField(Qp(5), [1, 0, 1])
   Traceback (most recent call last):
   ...
   padic_tools.core.DomainError: the defining polynomial must be irreducible modulo p; split and ramified presentations are unsupported
   >>> UnramifiedQuadraticField(Qp(2), [1, 0, 1])
   Traceback (most recent call last):
   ...
   padic_tools.core.DomainError: the defining polynomial must be irreducible modulo p; split and ramified presentations are unsupported

For an unramified quadratic extension of Q2 use a different polynomial:

.. doctest::

   >>> L2 = UnramifiedQuadraticField(Qp(2, 12), [1, 1, 1], name="w")
   >>> w = L2.generator
   >>> (w*w+w+1).is_zero()
   True
   >>> L2.residue_field().order
   4

Defining data and valuation
---------------------------

Coefficients are exact rationals in ascending order, or a ``Polynomial``.
Write f(X)=X^2+B*X+C. The constructor requires monicity, degree two and
p-integral coefficients. It then checks that the reduction of f is
irreducible over Fp. For odd p the discriminant must be a nonsquare; for
p=2 the reduction must have no root at either 0 or 1. Testing the
discriminant uses modular exponentiation rather than enumerating Fp.

Let theta be the chosen root. The reduction gives the field basis
1, theta-bar of F_(p**2), and the ring of integers has integral basis
1, theta. After dividing a+b*theta by the smallest coefficient power of p,
at least one reduced coefficient is nonzero. Independence of the residue
basis prevents cancellation. Thus

.. math::

   v_L(p)=1,\qquad
   v_L(a+b\theta)=\min\{v_p(a),v_p(b)\},\qquad
   \mathcal O_L=\mathbb Z_p+\mathbb Z_p\theta.

The uniformizer remains p: the ramification index is one and the residue
degree is two. These are the conventions of `Sage's unramified extensions
<https://doc.sagemath.org/html/en/reference/padics/sage/rings/padics/unramified_extension_generic.html>`_.
They also agree with `Magma's precision normalization
<https://docs.magma-maths.org/LocalFields/PAdicRingsAndTheirExtensions/overview.html>`_
when the ramification index is one. The coordinate-minimum formula must not
be reused for a ramified power basis.

``L(a,b)`` and ``L((a,b))`` mean a+b*theta. Exact nonzero input has the base
field's default number of significant digits. Zero input gives a zero ball
at the base precision. Explicit ``absprec=N`` selects an absolute precision,
including a negative N. Qp coefficient balls and existing elements may be
converted, but the result cannot exceed their common known precision.
Conversions between different defining polynomials need a separately chosen
isomorphism; even isomorphic abstract residue fields do not choose one for us.

.. doctest::

   >>> y = L(Fraction(1, 3), 2)
   >>> (y.valuation(), y.precision, y.relative_precision)
   (-1, 11, 12)
   >>> y.absolute_value()
   Fraction(3, 1)
   >>> limited = L(Qp(3)(1, absprec=4), Qp(3)(2, absprec=7))
   >>> tuple(c.precision for c in limited.coefficients)
   (4, 4)
   >>> L(limited, absprec=5)
   Traceback (most recent call last):
   ...
   padic_tools.core.PrecisionError: cannot manufacture digits of an existing ball

Ideal balls and error propagation
---------------------------------

Every element represents the whole ideal ball x0+p^N O_L. Both coordinate
balls have absolute precision N. An input rectangle with unequal coordinate
precisions is enclosed by the common ideal ball; the additional coordinate
digits are discarded explicitly. The model tracks no correlations between
repeated uses of a ball.

.. _extension-precision-contract:

For x and y of absolute precisions N and M and valuation lower bounds v and
w (a zero ball uses its absolute precision as that lower bound):

.. list-table:: Guaranteed absolute precision
   :header-rows: 1
   :widths: 40 60

   * - Operation
     - Resulting precision
   * - x+y, x-y
     - min(N,M)
   * - x*y
     - min(N+w,M+v)
   * - inverse(x), when zero is excluded
     - N-2*v
   * - conjugate(x), Frobenius(x), trace(x)
     - N
   * - algebraic norm(x)
     - N+v; for a zero ball, 2*N
   * - x times an exact nonzero scalar q
     - N+v_p(q)
   * - x divided by an exact nonzero scalar q
     - N-v_p(q)

For multiplication, write the perturbed inputs as x+delta and y+epsilon.
The three error terms x*epsilon, y*delta and delta*epsilon have valuations
at least v+M, w+N and N+M. Since v is at most N and w at most M, the last
term does not weaken the bound. For inversion, the difference of inverses
is -delta/(x*(x+delta)); both denominator factors have valuation v.

Conjugation preserves the integral ideal and is an isometry. With the
nontrivial automorphism sigma, trace(delta) stays in p^N Zp. The norm error
is trace(sigma(x)*delta)+norm(delta), with bounds N+v and 2*N. This proves
the stated norm bound, including the zero-ball case. For a nonzero x,
v_p(norm(x))=2*v_L(x), since both conjugates have the same valuation.

.. doctest::

   >>> x = L(9, 18, absprec=10)
   >>> (x.valuation(), x.inverse().precision, x.norm().precision)
   (2, 6, 12)
   >>> ((x*9).precision, (x/9).precision, x.trace().precision)
   (12, 8, 10)
   >>> z = L.O(-2)
   >>> z.norm().precision, z.contains(Fraction(1, 9), 0)
   (-4, True)
   >>> z.inverse()
   Traceback (most recent call last):
   ...
   padic_tools.core.PrecisionError: cannot invert an extension ball containing zero

``norm()`` here means the algebraic norm to Qp. Use ``absolute_value()``
for p**(-v); a zero ball has no known nonzero valuation or absolute value.
``contains`` accepts exact rational coordinates. ``is_congruent`` compares
balls at their common precision, while equality compares identical balls.
An equality of overlapping balls is not an equality of underlying exact
numbers. ``representative()`` returns one rational coordinate pair; it does
not select a complex image.

Residues and Frobenius
----------------------

``residue()`` requires integral input and at least one known absolute digit.
It returns an exact element of Fp[t]/(f-bar), with arithmetic, inverses,
integer powers, trace, norm and Frobenius. Equality there is exact.

.. doctest::

   >>> R = L.residue_field()
   >>> R.order
   9
   >>> r = L(4, 5).residue()
   >>> r.coefficients
   (1, 2)
   >>> r*r.inverse() == R(1)
   True
   >>> r**3 == r.frobenius() == r.conjugate()
   True
   >>> r**9 == r
   True
   >>> L.lift_residue(r).residue() == r
   True

``lift_residue`` chooses integer coefficient representatives of an exact
residue; it is a section of reduction, not a ring map to characteristic zero
and not a Teichmueller lift. Extra digits specify that chosen integer lift.
The lift of Frobenius on this unramified degree-two local field is the
nontrivial automorphism. In the selected basis,

.. math::

   \sigma(\theta)=-B-\theta,\qquad
   \operatorname{Tr}(a+b\theta)=2a-Bb,\qquad
   \operatorname{N}(a+b\theta)=a^2-Bab+Cb^2.

Residue Frobenius is the p-power map. Local-field Frobenius is
``conjugate()``, not the p-th power map in characteristic zero:

.. doctest::

   >>> x = L(1, 1)
   >>> x.frobenius().residue() == x.residue()**3
   True
   >>> x.frobenius().frobenius() == x
   True
   >>> x.frobenius().is_congruent(x**3)
   False

Relation to the bridge and remaining scope
------------------------------------------

``QuadraticAlgebra`` continues to describe a formal coefficient algebra,
including split and ramified cases. The new field constructor certifies
the unramified case and supplies its valuation and ideal precision. A common
exact Gaussian source can be evaluated in Q3(i) by explicitly choosing
theta=i:

.. doctest::

   >>> source = Gaussian(Fraction(1, 3), 2)
   >>> image = L(source.a, source.b)
   >>> image.norm().contains(source.norm())
   True
   >>> square = source*source
   >>> (image*image).contains(square.a, square.b)
   True

The existing ``GaussianBridge.padic_algebra`` and quaternion bridge still
return their documented coefficient algebras. The separate :doc:`ramified`
engine now includes Q2(i) in a certified uniformizer basis. The
:doc:`local_bridges` adapters select actual local fields with an explicit
generator image. Other extension degrees, general isomorphisms and extension
analytic functions remain later steps. Existing Qp analytic functions do
not acquire a larger domain merely by introducing this field class.
``Polynomial`` evaluation already works by arithmetic in the chosen field.

Validation includes exhaustive small residue fields, SymPy polynomial
quotient inverses and multiplication-matrix determinants/traces, perturbed
input-ball enclosure tests at positive and negative precision, and an
optional independent PARI/GP polymod oracle for p=2,3,5,7. The latter uses
PARI's Qp arithmetic and polynomial quotient inversion; it does not claim a
comparison against a general ramified local-field backend.
