"""Independent exact oracles, input-ball perturbations and conditioning."""
from dataclasses import FrozenInstanceError
from fractions import Fraction as F
from itertools import combinations, product
import random

import pytest
import sympy as sp

from padic_tools import (
    CliffordAlgebra, DomainError, EisensteinQuadraticField, Gaussian,
    GaussianLocalPlace, Hamilton, Matrix, PrecisionError, Q2i, Qp,
    QuadraticAlgebra, QuadraticElement, UnramifiedQuadraticField, exact, valuation,
)


t = sp.Symbol("t")


def expression(q):
    if isinstance(q, tuple):
        return sp.Rational(q[0])+sp.Rational(q[1])*t
    if hasattr(q, "parent") and q.parent == Gaussian:
        return sp.Rational(q.a)+sp.I*sp.Rational(q.b)
    return sp.Rational(q)


def reduced(q, f):
    numerator, denominator = sp.cancel(q).as_numer_denom()
    q = sp.rem(numerator*sp.invert(denominator, f, t), f, t)
    return F(q.coeff(t, 0)), F(q.coeff(t, 1))


def quotient_inverse(M, f):
    return [[reduced(c, f) for c in row] for row in M.inv().tolist()]


def gaussian(q):
    a, b = sp.expand(q).as_real_imag()
    return Gaussian(F(a), F(b))


@pytest.mark.parametrize("base", [exact, Gaussian])
@pytest.mark.parametrize("n", [1, 2, 3, 4])
def test_exact_matrices_against_sympy(base, n):
    rng = random.Random(129+n)
    for _ in range(3):
        rows = [[F(rng.randint(-5, 5), rng.randint(1, 4)) for _ in range(n)] for _ in range(n)]
        if base == Gaussian:
            rows = [[Gaussian(c, rng.randint(-2, 2)) for c in row] for row in rows]
        M = sp.Matrix([[expression(c) for c in row] for row in rows])
        A = Matrix(base, rows)
        converter = F if base is exact else gaussian
        assert A.determinant() == converter(M.det())
        assert A.trace() == converter(sp.trace(M))
        assert (A@A).rows == tuple(tuple(converter(c) for c in row) for row in (M*M).tolist())
        if M.det() == 0:
            with pytest.raises(ZeroDivisionError):
                A.inverse()
            continue
        assert A.inverse().rows == tuple(tuple(converter(c) for c in row) for row in M.inv().tolist())
        rhs = tuple(range(1, n+1))
        reference = M.inv()*sp.Matrix(rhs)
        info = A.solve_with_info(rhs)
        assert info.solution == tuple(converter(c) for c in reference)
        assert info.condition_exponent is info.base_condition_exponent is None
        assert info.pivot_valuations == () and info.precisions == (None,)*n
        assert A@A.inverse() == Matrix.identity(base, n)


@pytest.mark.parametrize("p", [2, 3, 5, 7])
@pytest.mark.parametrize("n", [1, 2, 3, 4])
def test_qp_matrices_and_complete_pivot_invariants_against_minors(p, n):
    rng = random.Random(200+p+n)
    for _ in range(3):
        rows = [[F(rng.randint(-6, 6), p**rng.randint(0, 2)) for _ in range(n)] for _ in range(n)]
        M = sp.Matrix(rows)
        A = Matrix(Qp(p, 20), rows, absprec=20)
        assert A.determinant().contains(F(M.det()))
        if M.det() == 0:
            with pytest.raises(PrecisionError):
                A.inverse()
            continue
        inverse = A.inverse()
        assert inverse.contains([[F(c) for c in row] for row in M.inv().tolist()])
        rhs = tuple(F(i, p) for i in range(1, n+1))
        info = A.solve_with_info(rhs)
        assert all(x.contains(F(y)) for x, y in zip(info.solution, M.inv()*sp.Matrix(rhs)))
        # Determinantal divisors provide an independent Smith-valuation oracle.
        previous = 0
        smith = []
        for size in range(1, n+1):
            divisors = [valuation(F(M.extract(rr, cc).det()), p)
                        for rr in combinations(range(n), size) for cc in combinations(range(n), size)]
            value = min(divisors)
            smith.append(value-previous)
            previous = value
        assert tuple(smith) == info.pivot_valuations
        assert info.condition_exponent == smith[-1]-smith[0]
        assert info.base_condition_exponent == info.condition_exponent
        assert sorted(info.row_permutation) == sorted(info.column_permutation) == list(range(n))


@pytest.mark.parametrize("p,f,weights", [
    (2, (1, 1, 1), (0, 0)), (3, (1, 0, 1), (0, 0)),
    (2, (2, -2, 1), (0, 1)), (3, (-3, 0, 1), (0, 1)),
])
@pytest.mark.parametrize("n", [2, 3])
def test_quadratic_solutions_and_conditions_against_sympy_quotients(p, f, weights, n):
    base = (UnramifiedQuadraticField if weights == (0, 0) else EisensteinQuadraticField)(Qp(p, 18), f)
    e = base.ramification_index
    modulus = sum(c*t**i for i, c in enumerate(f))
    rng = random.Random(p*200+n+sum(f))
    rows = [[(F(rng.randint(-4, 4), p), F(rng.randint(-3, 3), p)) for _ in range(n)] for _ in range(n)]
    M = sp.Matrix([[expression(c) for c in row] for row in rows])
    A = Matrix(base, rows, absprec=23)
    assert A.determinant().contains(reduced(M.det(), modulus))
    inverse = A.inverse()
    assert inverse.contains(quotient_inverse(M, modulus))
    rhs = tuple(range(1, n+1))
    info = A.solve_with_info(rhs)
    references = [reduced(c, modulus) for c in M.inv()*sp.Matrix(rhs)]
    assert all(x.contains(y) for x, y in zip(info.solution, references))
    previous, smith = 0, []
    for size in range(1, n+1):
        minors = [reduced(M.extract(rr, cc).det(), modulus)
                  for rr in combinations(range(n), size) for cc in combinations(range(n), size)]
        v = min(min(e*valuation(c, p)+w for c, w in zip(minor, weights)) for minor in minors)
        smith.append(v-previous)
        previous = v
    assert info.pivot_valuations == tuple(smith)
    assert info.base_condition_exponent == F(smith[-1]-smith[0], e)


@pytest.mark.parametrize("p,f", [(3, None), (3, (1, 0, 1)), (2, (2, -2, 1))])
@pytest.mark.parametrize("n", [-1, 1, 3])
def test_solutions_enclose_perturbations_of_entire_input_balls(p, f, n):
    if f is None:
        K, pi, modulus = Qp(p, 12), sp.Integer(p), None
        convert = lambda q: F(q)
    else:
        K = (UnramifiedQuadraticField if p == 3 else EisensteinQuadraticField)(Qp(p, 12), f)
        pi = sp.Integer(p) if p == 3 else t
        modulus = sum(c*t**i for i, c in enumerate(f))
        convert = lambda q: reduced(q, modulus)
    scale = pi**(n-2)
    center = sp.Matrix([[scale, 2*scale], [0, scale]])
    rhs = sp.Matrix([scale, -scale])
    A = Matrix(K, [[convert(c) for c in row] for row in center.tolist()], absprec=n)
    b = [K(convert(c), absprec=n+1) for c in rhs]
    computed = A.solve(b)
    inverse = A.inverse()
    determinant = A.determinant()
    for signs in product((-1, 1), repeat=4):
        perturbed = center+pi**n*sp.Matrix(2, 2, signs)
        assert A.contains([[convert(c) for c in row] for row in perturbed.tolist()])
        reference_inverse = perturbed.inv()
        assert inverse.contains([[convert(c) for c in row] for row in reference_inverse.tolist()])
        assert determinant.contains(convert(perturbed.det()))
        for db in ((1, -1), (-1, 1)):
            actual_rhs = rhs+pi**(n+1)*sp.Matrix(db)
            reference = reference_inverse*actual_rhs
            assert all(x.contains(convert(y)) for x, y in zip(computed, reference))


@pytest.mark.parametrize("p", [2, 3, 5])
def test_uncertain_zero_entries_cannot_be_discarded(p):
    K = Qp(p, 20)
    A = Matrix(K, [[K(1), K.O(2)], [K.O(2), K(1)]])
    inverse = A.inverse()
    reference = sp.Matrix([[1, p*p], [p*p, 1]]).inv()
    assert inverse.contains([[F(c) for c in row] for row in reference.tolist()])
    assert inverse[0, 0].precision == 4
    assert inverse[0, 1].precision == inverse[1, 0].precision == 2
    assert A.determinant().precision == 4


@pytest.mark.parametrize("base", [Qp(3, 3), Q2i(3), UnramifiedQuadraticField(Qp(3, 3), [1, 0, 1])])
def test_high_precision_outputs_do_not_inherit_parent_default(base):
    A = Matrix(base, [[1, 2], [3, 7]], absprec=19)
    assert min(c.precision for row in A.inverse() for c in row) >= 19
    assert (A@A).is_congruent(A**2, 19)
    assert all(c.precision == 19 for row in (A*0) for c in row)
    zeros = A.solve([0, 0])
    assert all(c.is_zero() and c.precision >= 19 for c in zeros)
    for op in (lambda: A.truncate(20), lambda: Matrix(base, A.rows, absprec=20),
               lambda: A.is_congruent(A, 20)):
        with pytest.raises(PrecisionError):
            op()
    assert A.truncate(7).precisions == ((7, 7), (7, 7))
    mixed = Matrix(base, [[base(2, absprec=19), base(2, absprec=19)],
                          [base(3, absprec=19), base(7, absprec=3)]])
    with pytest.raises(PrecisionError):
        A.is_congruent(mixed, 17)  # A first-entry mismatch cannot hide missing later digits.
    scale = 2 if base.p == 2 else 3
    assert (A*F(1, scale))[0, 0].precision == 19-(2 if base.p == 2 else 1)


@pytest.mark.parametrize("base", [Qp(3, 12), Q2i(12)])
def test_condition_exponent_and_absolute_inverse_loss(base):
    pi = base(3, absprec=32) if isinstance(base, type(Qp(3))) else base(0, 1, absprec=32)
    A = Matrix(base, [[pi**2, 0], [0, pi**5]], absprec=21)
    info = A.solve_with_info([1, 1])
    assert info.pivot_valuations == (2, 5) and info.condition_exponent == 3
    assert info.base_condition_exponent == F(3, 1 if base.p == 3 else 2)
    inverse = A.inverse()
    assert inverse[0, 0].precision == 17 and inverse[1, 1].precision == 11
    scaled = (A*F(1, 3 if base.p == 3 else 2)).solve_with_info([1, 1])
    assert scaled.condition_exponent == 3


def test_complete_column_pivoting_restores_original_unknown_order():
    K = Qp(3, 15)
    A = Matrix(K, [[9, 1], [3, 0]], absprec=15)
    info = A.solve_with_info([11, 3])
    assert info.column_permutation == (1, 0)
    assert all(x.contains(y) for x, y in zip(info.solution, (1, 2)))
    rhs = Matrix(K, [[11, -1], [3, 0]], absprec=15)
    solution = A.solve(rhs)
    assert isinstance(solution, Matrix) and solution.contains([[1, 0], [2, -1]])


@pytest.mark.parametrize("p", [2, 3, 5])
def test_matrix_base_change_commutes_with_products_and_solutions(p):
    place = GaussianLocalPlace(Qp(p, 20), root_sign=-1)
    A = Matrix(Gaussian, [[Gaussian(1, 1), 2], [3, Gaussian(2, -1)]])
    B = Matrix(Gaussian, [[1, Gaussian(1, -2)], [0, 1]])
    mapper = lambda c: place(c, absprec=19)
    local = A.map_coefficients(place.field, mapper)
    assert (local@B.map_coefficients(place.field, mapper)).is_congruent((A@B).map_coefficients(place.field, mapper))
    assert local.inverse().is_congruent(A.inverse().map_coefficients(place.field, mapper))
    rhs = Matrix(Gaussian, [[1], [Gaussian(2, 1)]])
    assert local.solve(rhs.map_coefficients(place.field, mapper)).is_congruent(A.solve(rhs).map_coefficients(place.field, mapper))


def test_shapes_exact_errors_immutability_and_domain_checks():
    A = Matrix(exact, [[1, 2, 3], [4, 5, 6]])
    assert A.shape == (2, 3) and A.T.shape == (3, 2)
    assert A.T.T == A and A[1, 2] == 6 and A[0] == (1, 2, 3)
    assert (A+A).rows == (A*2).rows and (A-A).rows == (A*0).rows
    assert (2*A).rows == (A*2).rows and (-A+A).rows == (A*0).rows
    with pytest.raises(FrozenInstanceError):
        A.rows = ()
    for op in (A.determinant, A.trace, A.inverse, lambda: A.solve([1, 2]), lambda: A**2):
        with pytest.raises(ValueError):
            op()
    for rows in ([], [[]], [[1], [1, 2]]):
        with pytest.raises(ValueError):
            Matrix(exact, rows)
    with pytest.raises(DomainError):
        Matrix(exact, [[1]*65])
    for base in (object(), Hamilton, CliffordAlgebra(), QuadraticAlgebra(Qp(3))):
        with pytest.raises(TypeError):
            Matrix(base, [[1]])
    for value in (1.0, 1j, True):
        with pytest.raises(TypeError):
            Matrix(exact, [[value]])
    with pytest.raises(TypeError):
        Matrix(Gaussian, [[QuadraticElement(Gaussian, Qp(3)(1), 0)]])
    with pytest.raises(TypeError):
        Matrix(exact, [[1]], absprec=5)
    for n in (True, 1.5, "2"):
        with pytest.raises(TypeError):
            Matrix(Qp(3), [[1]], absprec=n)
    for size in (0, -1, True, 1.5, 65):
        with pytest.raises(ValueError):
            Matrix.identity(exact, size)
    square = Matrix(exact, [[1, 2], [3, 7]])
    assert square**0 == Matrix.identity(exact, 2)
    assert square**(-1) == square.inverse()
    for op in (lambda: square+Matrix(Qp(3), [[1, 2], [3, 7]]),
               lambda: square@2, lambda: square* square, lambda: square**True,
               lambda: square.truncate(3), lambda: square.is_congruent(square, 3)):
        with pytest.raises(TypeError):
            op()
    for op in (lambda: square+A, lambda: A@square, lambda: square.solve([1]),
               lambda: square.solve(Matrix(exact, [[1]])), lambda: square.contains([[1]])):
        with pytest.raises(ValueError):
            op()


@pytest.mark.parametrize("p", [2, 3, 5])
def test_insufficient_precision_is_not_an_exact_singularity_claim(p):
    K = Qp(p, 12)
    for rows in ([[0, 0], [0, 0]], [[1, 1], [1, 1]],
                 [[K(p**2), K.O(0)], [K.O(0), K(p**3)]]):
        A = Matrix(K, rows)
        # Determinants remain available even if an inverse cannot be certified.
        A.determinant()
        with pytest.raises(PrecisionError):
            A.inverse()
    singular = Matrix(exact, [[1, 1], [1, 1]])
    assert singular.determinant() == 0
    with pytest.raises(ZeroDivisionError):
        singular.solve([1, 1])
    with pytest.raises(TypeError):
        Matrix(K, [[1]]).contains([[K(1)]])
    with pytest.raises(TypeError):
        Matrix(K, [[1, 1]]).contains([[2, K(1)]])
    with pytest.raises(TypeError):
        Matrix(K, [[1]]).solve([Qp(7)(1)])
