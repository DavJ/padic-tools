Matrices and precision-aware linear systems
===========================================

``Matrix(base, rows)`` supplies immutable dense matrices over exact Q,
exact Q[i], Qp, and either certified quadratic local-field engine. The
native implementation uses the standard library. SymPy and PARI/GP are
independent test oracles, rather than runtime dependencies.

Exact sources and basic operations
----------------------------------

Use ``exact`` for rational entries and ``Gaussian`` for rational Gaussian
entries. Matrix multiplication is ``@``; ``*`` scales by a scalar. Addition,
subtraction, integer powers, transpose, trace and determinant are available.
``T`` transposes without conjugating coefficients. Matrices are nonempty
and rectangular; each dimension is limited to 64 in the native backend.
Trace, determinant, inverse and solves require a square matrix.

.. doctest::

   >>> A = Matrix(exact, [[1, 2], [3, 7]])
   >>> A.shape, A.determinant(), A.trace()
   ((2, 2), Fraction(1, 1), Fraction(8, 1))
   >>> A.solve([1, 0])
   (Fraction(7, 1), Fraction(-3, 1))
   >>> A @ A.inverse() == Matrix.identity(exact, 2)
   True
   >>> A.T.T == A and A**2 == A@A and A*2 == A+A
   True
   >>> B = Matrix(exact, [[1, 2], [0, 1]])
   >>> X = A.solve(B)  # A Matrix right-hand side returns a Matrix.
   >>> A@X == B
   True
   >>> G = Matrix(Gaussian, [[Gaussian(1, 1), 2], [3, Gaussian(2, -1)]])
   >>> G @ G.inverse() == Matrix.identity(Gaussian, 2)
   True

A vector right-hand side returns a tuple. Multiple right-hand sides must
be passed as a ``Matrix``; nested sequences otherwise describe the scalar
vector entries, including coordinate pairs in a quadratic field. Bases
must agree for matrix arithmetic. Use ``map_coefficients`` for an explicit
base change, or construct a new matrix over the intended base.

Local matrices are whole input families
---------------------------------------

Each entry in a local matrix is a scalar ball. The matrix represents their
Cartesian product, allowing entries to vary independently. Repeated use of
one ball does not record a correlation. In particular, a zero-containing
entry is an uncertainty, rather than an exact structural zero.

``absprec=N`` sets every entry's absolute precision in uniformizer units.
Without it, the scalar field's normal construction rules apply. Existing
balls retain their precision, and cannot acquire digits through conversion
or ``truncate``. In an Eisenstein quadratic field the unit of precision is
pi, with v_pi(p)=2. In Qp and unramified fields it is p.

.. doctest::

   >>> K = Qp(3, 3)
   >>> local = Matrix(K, [[1, 2], [3, 7]], absprec=19)
   >>> local.precisions
   ((19, 19), (19, 19))
   >>> result = local.solve([1, 0])
   >>> all(x.contains(y) for x, y in zip(result, (7, -3)))
   True
   >>> local.inverse().contains(A.inverse())
   True
   >>> local.truncate(7).precisions
   ((7, 7), (7, 7))
   >>> local.solve([0, 0])[0].precision >= 19
   True

Exact rational vector right-hand sides stay exact inside elimination.
An inverse uses an internally exact identity, so it is not capped by a
field default below the matrix's known precision. By contrast,
``Matrix.identity(K, n)`` constructs finite balls for all its entries;
choose ``absprec`` explicitly when using it as input. ``M**0`` returns
that default-precision identity. Positive powers avoid introducing one.

For a nonrational local right-hand side, construct scalar balls or a
right-hand-side matrix at the needed precision. ``contains`` checks exact
entries, not p-adic centres extracted from approximate input. Use
``is_congruent(other, precision=N)`` for comparisons of two ball matrices;
it raises if any entry lacks the requested digits.

Certified complete pivoting
---------------------------

The square nonsingular solver chooses a nonzero pivot of minimum valuation
in the entire remaining block, and swaps both rows and columns. It also
checks the lower bounds of zero-containing balls. A ball with valuation
bound below the chosen pivot could hide a larger p-adic entry; the solver
then raises ``PrecisionError`` instead of guessing its valuation.

For a successful solve, every matrix in the represented input family is
invertible, and every corresponding solution lies in the returned entry
balls. The proof follows the scalar enclosure rules through elimination:
each chosen pivot excludes zero, so division is valid for every represented
input. Row operations and back substitution enclose their exact results.
The eliminated entry is structurally zero by the row-operation identity;
other zero-containing entries are still processed. The final column
permutation is undone to restore the original order of unknowns.

.. doctest::

   >>> P = Matrix(Qp(3, 15), [[9, 1], [3, 0]], absprec=15)
   >>> info = P.solve_with_info([11, 3])
   >>> info.column_permutation
   (1, 0)
   >>> all(x.contains(y) for x, y in zip(info.solution, (1, 2)))
   True

If no pivot can be certified, local matrices raise ``PrecisionError``.
An exactly singular matrix over Q or Q[i] raises ``ZeroDivisionError``.
A singular rational centre alone does not prove that a local ball matrix
is singular. Conversely, an invertible centre alone does not certify all
matrices in its input family. Reevaluate an exact source at higher
precision when the available balls are insufficient.

Condition exponents and output precision
----------------------------------------

Let s_1,...,s_n be the successful pivot valuations. They are nondecreasing
elementary-divisor valuations, in uniformizer units. To see this, the
minimum-valuation test makes all row-elimination coefficients integral.
The remaining upper row can likewise be cleared by integral column
operations. These operations and their inverses are unimodular over the
valuation ring; repeating on the remaining block gives a diagonal form
with the displayed ordered valuations.

Unimodular transformations preserve the vector sup norm. If e is the
ramification index, the matrix sup-norm condition number therefore is

.. math::

   \kappa_\infty(A)=\lVert A\rVert_\infty\lVert A^{-1}\rVert_\infty
      =p^{(s_n-s_1)/e}.

``LinearSolveResult.condition_exponent`` is s_n-s_1;
``base_condition_exponent`` is the exact rational (s_n-s_1)/e. These
diagnostics are certified for the entire input family when elimination
succeeds. They are ``None`` for an exact base. The permutation tuples list
original indices in elimination order, while ``solution`` is already
returned in original order.

The condition exponent describes relative sensitivity, not one universal
number to subtract from every absolute precision. Input scale, right-hand
side precision and individual entries also matter. For example, uncertainty
in a diagonal a+O(pi^N) gives inverse precision N-2*v(a). Uncertain
off-diagonal entries also affect a solution's precision.

.. doctest::

   >>> D = Matrix(Qp(3, 12), [[9, 0], [0, 243]], absprec=21)
   >>> report = D.solve_with_info([1, 1])
   >>> report.pivot_valuations, report.condition_exponent, report.base_condition_exponent
   ((2, 5), 3, Fraction(3, 1))
   >>> D.inverse().precisions
   ((17, 14), (14, 11))
   >>> report.precisions
   (14, 11)
   >>> L2 = Q2i(12)
   >>> pi = L2(0, 1, absprec=32)
   >>> D2 = Matrix(L2, [[pi**2, 0], [0, pi**5]], absprec=21)
   >>> D2.solve_with_info([1, 1]).base_condition_exponent
   Fraction(3, 2)

The native model propagates separate balls. It gives safe enclosures and
may lose precision through repeated dependencies; it does not implement
optimal lattice precision. See `Caruso, Roe and Vaccon, p-Adic Stability
in Linear Algebra <https://arxiv.org/abs/1506.05644>`_ for that broader
precision framework, and `Kulkarni, section 2
<https://arxiv.org/abs/1907.03740>`_ for p-adic conditioning and normal forms.
The complete-pivot certificate above applies the valuation-ring argument
also to the library's certified quadratic fields.

Determinants without division
-----------------------------

``determinant()`` uses the division-free Berkowitz recurrence, with O(n^4)
scalar arithmetic. Unlike solving, it needs neither a pivot nor a nonzero
test, and remains available for singular or unresolved matrices.

For a block matrix with first row/column a, r, c and trailing block B,
form the Toeplitz column ``1, -a, -r*c, -r*B*c, ...``. Multiplication of
its lower-triangular Toeplitz matrix by the descending coefficients of
det(tI-B) gives those of det(tI-A). Starting with the empty-block
polynomial 1 yields the determinant as (-1)^n times the constant term.
Only additions and multiplications of ball enclosures are used. See the
`Sage matrix reference
<https://doc.sagemath.org/html/en/reference/matrices/sage/matrix/matrix2.html>`_
for the division-free setting.

.. doctest::

   >>> singular_source = Matrix(exact, [[1, 1], [1, 1]])
   >>> singular_source.determinant()
   Fraction(0, 1)
   >>> unresolved = singular_source.map_coefficients(Qp(3, 10))
   >>> unresolved.determinant().is_zero()
   True

The last result is a determinant ball containing zero, not an exact
singularity certificate. Rank, kernels, rectangular solves, eigenvalues,
noncommutative coefficient solvers and sparse/accelerated backends remain
separate extensions. The solver reports invariant valuations, rather than
returning a full Smith decomposition with transformation matrices.

Compatible Gaussian places
--------------------------

Exact source matrices can be mapped through any :doc:`local_bridges`
Gaussian embedding. Matrix products, determinants and defined inverses
commute with that scalar map, up to the propagated local ball precision.

.. doctest::

   >>> place = GaussianLocalPlace(Qp(2, 12), root_sign=-1)
   >>> mapper = lambda c: place(c, absprec=23)
   >>> local_G = G.map_coefficients(place.field, mapper)
   >>> local_G.inverse().is_congruent(G.inverse().map_coefficients(place.field, mapper))
   True
   >>> local_G.determinant().is_congruent(mapper(G.determinant()))
   True
   >>> (local_G@local_G).is_congruent((G@G).map_coefficients(place.field, mapper))
   True

Validation and a runnable example
---------------------------------

Tests compare exact inverses and products with SymPy, local determinants
and solves with independent rational/quadratic quotient calculations, and
pivot valuations with valuations of all determinantal minors in small
matrices. Whole-input perturbations exercise negative and odd uniformizer
precisions. PARI/GP checks Qp and both Gaussian completion types. The
zero-ball regression changes unresolved off-diagonal entries and verifies
their effect on inverse diagonal precision.

Run ``python examples/linear_systems.py`` for exact, conditioned Q3 and
ramified Q2(i) examples without optional dependencies.
