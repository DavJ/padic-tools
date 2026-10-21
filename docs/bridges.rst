Complex, p-adic, quaternion and Clifford bridges
================================================

The implemented bridge preserves a common exact source and records which
completion or coefficient algebra receives it. This chapter distinguishes
algebraic identifications from analytic comparison problems.

The diagonal rational bridge
----------------------------

An exact q in Q has an archimedean image and a p-adic image:

.. math::

   \mathbb Q\longrightarrow\mathbb C\times\mathbb Q_p,
   \qquad q\longmapsto(\iota_\infty(q),\iota_p(q)).

Both maps preserve sums and products. The complex image is rounded using an
isolated mpmath context, and the p-adic image is a ball with known precision.
The exact rational remains the source of truth and can be reevaluated.

.. doctest::

   >>> source = RationalBridge(Fraction(2, 7))
   >>> source.padic(Qp(5)).contains(Fraction(2, 7))
   True
   >>> source.value
   Fraction(2, 7)

There is no continuous field map Qp -> C fixing Q. Indeed p^n tends to zero
in Qp; its image p^n does not tend to zero in C. This short obstruction is
why a finite p-adic digit string is not converted to a purported complex
approximation. ``representative()`` explicitly returns one rational centre
of a ball, which may be chosen in many inequivalent ways.

General number-field embeddings into different completions require choices
of places; compare `Sage's number-field embedding interface
<https://doc.sagemath.org/html/en/reference/number_fields/sage/rings/number_field/number_field_morphisms.html>`_.
Our native implementation provides Q and Q[i], rather than an unimplemented
claim of arbitrary number fields.

Q[i] and chosen local places
----------------------------

``Gaussian(a,b)`` stores the exact algebraic element a+b*i with i^2=-1.
``GaussianBridge`` supplies:

.. list-table::
   :header-rows: 1
   :widths: 27 73

   * - Method
     - Meaning
   * - ``complex(dps=...)``
     - The embedding i -> the usual complex imaginary unit, numerically rounded
   * - ``padic_algebra(K)``
     - Coefficient base change Q[i] -> Qp[t]/(t^2+1), with coefficient ball errors
   * - ``split_padic(K, root_sign=...)``
     - A chosen place i -> +sqrt(-1) or -sqrt(-1), when the root lies in Qp

For p congruent to 1 modulo 4, -1 is a square in Qp. The two split place
choices are explicit, and the quadratic quotient is Qp x Qp, not a field.
For p congruent to 3 modulo 4 it is an unramified quadratic field. At p=2
it is a ramified quadratic field. In every case the coefficient algebra is
implemented; a full extension-field valuation/precision engine is not.

.. doctest::

   >>> g = GaussianBridge(Fraction(1, 3), 2)
   >>> g.padic_algebra(Qp(3)).parent.d
   Fraction(-1, 1)
   >>> i5 = GaussianBridge(0, 1).split_padic(Qp(5), root_sign=1)
   >>> j5 = GaussianBridge(0, 1).split_padic(Qp(5), root_sign=-1)
   >>> (i5*i5).is_congruent(-1) and i5.is_congruent(-j5)
   True

Inversion uses the quadratic norm a^2+b^2. Zero divisors at split primes
and balls whose norm is unresolved at zero are rejected. The coefficient
name ``i`` is an algebra generator and does not make this quotient the
completed algebraic closure Cp.

Quaternion and biquaternion base change
---------------------------------------

``QuaternionAlgebra(base, a, b)`` imposes i^2=a, j^2=b and ij=-ji=k.
``Hamilton`` uses (-1,-1) over Q. ``Biquaternions`` uses the same parameters
over exact Q[i], giving eight rational components.

.. math::

   H_{\mathbb Q}\otimes_{\mathbb Q}\mathbb Q_p,
   \qquad H_{\mathbb Q}\otimes_{\mathbb Q}\mathbb Q[i].

The multiplication table is unchanged under scalar base change. Reduced
norm and conjugation are implemented, with inversion only when the norm
is invertible. A p-adic Hamilton algebra is not assumed to be a division
algebra: with these parameters it is split for odd p and division at p=2.
For background on quaternion algebras over different scalar fields, see
`Sage's quaternion algebra reference
<https://doc.sagemath.org/html/en/reference/algebras/sage/algebras/quatalg/quaternion_algebra.html>`_.

.. doctest::

   >>> q = Hamilton(1, 2, 3, 4)
   >>> q.norm()
   Fraction(30, 1)
   >>> q*q.inverse() == Hamilton(1)
   True
   >>> b = Biquaternions(Gaussian(1, 2), Gaussian(3, -1), 0, 1)
   >>> b5 = QuaternionBridge(b).padic(Qp(5))
   >>> (b5*b5).is_congruent(QuaternionBridge(b*b).padic(Qp(5)))
   True

``QuaternionBridge.complex`` returns a quaternion with four complex
coefficients, not a single complex number. Its ``padic`` method returns
four Qp coefficients or four coefficients in Qp[t]/(t^2+1), according to
the exact source. Arbitrary real/complex quaternion samples have no implicit
exact p-adic image.

The geometric-algebra identification
------------------------------------

For the diagonal Clifford algebra Cl(3,0), let e1^2=e2^2=e3^2=1 and distinct
basis vectors anticommute. The even part has basis 1,e2e3,e3e1,e1e2.
Our explicit convention is

.. math::

   i\mapsto e_2e_3,\qquad
   j\mapsto e_3e_1,\qquad
   k\mapsto -e_1e_2.

Every displayed bivector squares to -1, and
``(e2*e3)*(e3*e1)=-e1*e2``. Thus the quaternion relations hold and the map
is an isomorphism of four-dimensional algebras. It works over Q, Qp,
Q[i], and their stated coefficient base changes. Tensoring with a complex
coefficient algebra gives the biquaternion version. The library also
implements other diagonal Clifford metrics; see `the general Clifford
construction in Sage <https://doc.sagemath.org/html/en/reference/algebras/sage/algebras/clifford_algebra.html>`_.

.. doctest::

   >>> a, b = Hamilton(1,2,3,4), Hamilton(4,3,2,1)
   >>> (a*b).to_clifford() == a.to_clifford()*b.to_clifford()
   True
   >>> a.to_clifford().to_quaternion() == a
   True
   >>> C = CliffordAlgebra((1,1,1))
   >>> e1, e2, e3 = (C.e(k) for k in range(3))
   >>> (e1^e2) == e1*e2
   True
   >>> e1^e1 == C(0)
   True
   >>> C.pseudoscalar()*C.pseudoscalar() == C(-1)
   True

Vectors, bivectors, trivectors/pseudoscalars and grade projections are
available as algebraic objects. A differential operator, metric field,
connection, action or Hamiltonian needs additional explicit definitions;
the scalar bridge does not derive those physical structures.

What commutes, and what remains open
------------------------------------

.. list-table::
   :header-rows: 1
   :widths: 37 63

   * - Relationship
     - Status
   * - Rational algebra -> independent C and Qp images
     - Implemented; ring operations preserved, with stated approximation errors
   * - Q[i] -> chosen split place / quadratic quotient
     - Implemented, with place and split-field distinctions explicit
   * - Quaternion multiplication -> even Clifford multiplication
     - Explicit isomorphism, tested on basis elements and general coefficients
   * - Algebraic zeta special values -> KL interpolation
     - Implemented, including Euler factor and branch restriction
   * - Same finite q-series polynomial -> different completions
     - Its rational coefficients can be evaluated in either domain
   * - Analytic theta modular transformations across completions
     - No general theorem or API is claimed
   * - Complex zeta zeros -> p-adic zeta zeros
     - No correspondence follows from the implemented special-value relation

The useful next research bridge is therefore an explicit arithmetic family
with compatible coefficient data and local realizations, together with a
proved statement about the chosen invariant. Geometric algebra organizes
the noncommutative part; comparison of limits and zeros still requires its
own analytic argument.
