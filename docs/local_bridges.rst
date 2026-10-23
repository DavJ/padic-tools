Chosen places and local-field bridges
=====================================

``GaussianLocalPlace`` evaluates a common exact source in Q[i] in a
certified local field. It records the prime, the image of i, and the target
field. ``RationalBridge.local`` and ``QuaternionBridge.local`` use the same
target fields. These adapters complement the coefficient-algebra methods
in :doc:`bridges`.

Choosing the image of i
-----------------------

The keyword ``root_sign`` is required. For a sign s equal to +1 or -1,
the target and generator image are:

.. list-table::
   :header-rows: 1
   :widths: 25 25 30 20

   * - Prime
     - Target
     - Image of i
     - (e, f)
   * - p = 1 mod 4
     - Qp
     - s*r
     - (1, 1)
   * - p = 3 mod 4
     - Qp(t)
     - s*t
     - (1, 2)
   * - p = 2
     - Q2(pi)
     - s*(pi-1)
     - (2, 1)

Here e is the ramification index and f the residue degree over Qp.
The split root r satisfies r^2=-1 and lifts the smaller root modulo p.
The unramified generator t satisfies t^2+1=0, while pi^2-2*pi+2=0
defines the ramified uniformizer. ``place.field`` is respectively Qp,
``UnramifiedQuadraticField`` or ``EisensteinQuadraticField``.
In the nonsplit cases the two signs give conjugate embeddings into the
same presented completion, rather than two different primes of Q[i].
At split primes they choose the two distinct places. Explicit generator
images follow the convention of `Sage's number-field embeddings
<https://doc.sagemath.org/html/en/reference/number_fields/sage/rings/number_field/number_field_morphisms.html>`_.

.. doctest::

   >>> places = [GaussianLocalPlace(Qp(p, 12), root_sign=1) for p in (2, 3, 5)]
   >>> [(v.kind, v.ramification_index, v.residue_degree) for v in places]
   [('ramified', 2, 1), ('unramified', 1, 2), ('split', 1, 1)]
   >>> all((v.generator_image()**2).is_congruent(-1) for v in places)
   True
   >>> source = GaussianBridge(Fraction(1, 3), 2)
   >>> z2, z3, z5 = [source.local(v) for v in places]
   >>> z2.contains(Fraction(-5, 3), 2)  # a+b*i = (a-b)+b*pi
   True
   >>> z3.contains(Fraction(1, 3), 2)
   True
   >>> z5.field == Qp(5, 12)
   True

Every image is a whole local-field ball, with valuation, inverse and the
precision rules of its field engine. See :doc:`extensions` and
:doc:`ramified`. Python ``complex`` values and p-adic approximations are
not accepted as common exact sources.

Absolute precision and reevaluation
-----------------------------------

For all ``local`` methods and calls to a place, ``absprec=N`` specifies
absolute precision in powers of the target uniformizer: p^N at odd primes,
pi^N in Q2(i). The default is ``place.field.precision``: the base default
at odd primes, twice the base default at p=2. Negative N is allowed.
These adapters choose a fixed absolute output precision; direct native
field construction from exact scalars normally uses relative precision.

The exact source can be reevaluated above the field's default precision.
Truncating an existing ball cannot recover discarded digits. At split
primes, a+b*r needs r to precision at least N-v_p(b). The implementation
lifts the chosen root to ``max(1, N-v_p(b))`` before multiplying by the
exact b and adding a. At p=2 it first changes exact coordinates to
(a-s*b, s*b), then rounds in the uniformizer basis. This accounts for
nonintegral coefficients and cancellation before rounding.

.. doctest::

   >>> v5 = GaussianLocalPlace(Qp(5, 3), root_sign=1)
   >>> cancelling = GaussianBridge(-7, 1)
   >>> cancelling.local(v5, absprec=2).is_zero()
   True
   >>> fine = cancelling.local(v5, absprec=8)
   >>> (fine.precision, fine.valuation())
   (8, 2)
   >>> (fine*fine.inverse()).is_congruent(1)
   True
   >>> v2 = GaussianLocalPlace(Qp(2, 3), root_sign=-1)
   >>> image = GaussianBridge(Fraction(1, 64), Fraction(-1, 16)).local(v2, absprec=19)
   >>> image.precision, tuple(c.precision for c in image.coefficients)
   (19, (10, 9))
   >>> RationalBridge(Fraction(1, 8)).local(v2, absprec=7).contains(Fraction(1, 8))
   True

Conjugations and norms
----------------------

``place.conjugate()`` changes the image of the central Gaussian i. For
exact z=a+b*i, the two images multiply to the image of a^2+b^2:

.. math::

   \phi_s(z)\,\phi_{-s}(z)=a^2+b^2.

In nonsplit completions the opposite embedding agrees with local-field
conjugation; the algebraic ``norm()`` and ``trace()`` return Qp balls.
At a split place the target is Qp itself. There is no nontrivial Qp
automorphism swapping the roots within one embedding; both places are
needed for the Gaussian norm. ``Padic.norm()`` is the p-adic absolute
value, not this algebraic norm.

.. doctest::

   >>> g = Gaussian(1, 2)
   >>> v3 = GaussianLocalPlace(Qp(3, 12), root_sign=1)
   >>> v3(g).conjugate().is_congruent(v3.conjugate()(g))
   True
   >>> v3(g).norm().contains(5) and v3(g).trace().contains(2)
   True
   >>> (v5(g)*v5.conjugate()(g)).is_congruent(5)
   True

Quaternion ``conjugate()`` instead negates the three vector coefficients
and leaves their scalar Gaussian values unchanged. The reduced quaternion
norm lies in that scalar field. Applying a further field norm to Qp is a
separate operation.

Quaternions and geometric algebra
---------------------------------

Rational quaternion sources can target any certified field directly, or a
Gaussian place. Biquaternions require a ``GaussianLocalPlace`` to specify
how their central i is mapped. The parameters of ``QuaternionAlgebra``
are preserved. Reduced norms, multiplication and quaternion conjugation
therefore commute with scalar base change. Inverses commute wherever
defined and resolved by the input precision. See `Sage's quaternion algebras
<https://doc.sagemath.org/html/en/reference/algebras/sage/algebras/quatalg/quaternion_algebra.html>`_.

.. doctest::

   >>> R = EisensteinQuadraticField(Qp(3, 12), [-3, 0, 1])
   >>> QuaternionBridge(Hamilton(1, 2, 3, 4)).local(R).norm().contains(30)
   True
   >>> q = Biquaternions(Gaussian(1, 2), Gaussian(3, -1), 0, 1)
   >>> r = Biquaternions(2, Gaussian(-1, 3), 1, 0)
   >>> Q, Rq = QuaternionBridge(q).local(v2), QuaternionBridge(r).local(v2)
   >>> (Q*Rq).is_congruent(QuaternionBridge(q*r).local(v2))
   True
   >>> Q.norm().is_congruent(v2(q.norm()))
   True
   >>> Q.conjugate().is_congruent(QuaternionBridge(q.conjugate()).local(v2))
   True

The Hamilton-to-even-Clifford map in :doc:`bridges` commutes with these
scalar maps. ``Multivector.map_coefficients`` preserves the diagonal metric.
The central Gaussian i stays a scalar; the three anticommuting quaternion
units map to bivectors.

.. doctest::

   >>> C2 = CliffordAlgebra((1, 1, 1), v2.field)
   >>> mapped = q.to_clifford().map_coefficients(C2, v2)
   >>> mapped.is_congruent(Q.to_clifford())
   True
   >>> mapped.to_quaternion().is_congruent(Q)
   True
   >>> ((q*r).to_clifford().map_coefficients(C2, v2)).is_congruent(mapped*Rq.to_clifford())
   True

An absent sparse blade is exact zero. ``C2({})`` is that sparse zero,
whereas ``C2(0)`` stores a scalar zero ball at the field's default
precision. Stored zero balls retain uncertainty and are not pruned.
Comparisons with absent blades and conversion back to quaternion preserve
known coefficient precision, including requests above the default.

The adapters do not make every nonzero quaternion invertible: split
algebras have zero divisors, and a reduced norm ball containing zero cannot
be inverted. Extension-field analytic functions and comparisons of complex
and p-adic limits or zeta zeros require further mathematics.

Independent checks
-------------------

Tests cover both orientations, the three completion types, nonintegral
sources, negative/elevated precision, cancellation, norms, inverses and the
Clifford commuting square. A separate SymPy 2x2 matrix model checks Hamilton
products and reduced norms over Q[i]. Optional PARI/GP tests compare split
roots, quotient coordinates, inverse, trace and norm at higher precision;
CI installs GP for those checks.
