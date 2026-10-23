"""Chosen Gaussian places and scalar base changes, with independent oracles."""
from fractions import Fraction as F
import random

import pytest
import sympy as sp

from padic_tools import (
    Biquaternions, CliffordAlgebra, DomainError, EisensteinQuadraticField,
    Gaussian, GaussianBridge, GaussianLocalPlace, Hamilton, PrecisionError,
    Qp, QuadraticAlgebra, QuadraticElement, Quaternion, QuaternionAlgebra,
    QuaternionBridge, RationalBridge, UnramifiedQuadraticField,
)


@pytest.mark.parametrize("p,kind,e,f", [
    (2, "ramified", 2, 1), (3, "unramified", 1, 2),
    (5, "split", 1, 1), (7, "unramified", 1, 2), (13, "split", 1, 1),
])
@pytest.mark.parametrize("sign", [-1, 1])
def test_place_classification_generator_and_orientation(p, kind, e, f, sign):
    place = GaussianLocalPlace(Qp(p, 10), root_sign=sign)
    assert (place.kind, place.ramification_index, place.residue_degree) == (kind, e, f)
    assert place.field.precision == e*10
    i = place.generator_image(absprec=17)
    assert i.precision == 17 and (i*i).is_congruent(-1)
    assert i.is_congruent(-place.conjugate().generator_image(absprec=17))
    assert place.conjugate().conjugate() == place
    if kind == "split":
        small = min(r for r in range(p) if (r*r+1) % p == 0)
        assert i.residue(1) == sign*small % p
    elif kind == "unramified":
        assert i.contains(0, sign)
    else:
        assert i.contains(-sign, sign)
        assert place.field.generator.valuation() == 1


@pytest.mark.parametrize("p", [2, 3, 5, 7, 13])
@pytest.mark.parametrize("sign", [-1, 1])
def test_exact_gaussian_maps_preserve_operations_and_conjugate_norm(p, sign):
    place = GaussianLocalPlace(Qp(p, 16), root_sign=sign)
    rng = random.Random(p+100*sign)
    for _ in range(8):
        a, b = [Gaussian(F(rng.randint(-8, 8), p**rng.randint(0, 2)),
                         F(rng.randint(-8, 8), p**rng.randint(0, 2))) for _ in range(2)]
        A, B = place(a), place(b)
        assert (A+B).is_congruent(place(a+b))
        assert (A*B).is_congruent(place(a*b))
        if a.a or a.b:
            assert A.inverse().is_congruent(place(a.inverse()))
        opposite = place.conjugate()(a)
        assert opposite.is_congruent(place(a.conjugate()))
        assert (A*opposite).is_congruent(RationalBridge(a.norm()).local(place))
        if place.kind != "split":
            assert A.conjugate().is_congruent(opposite)
            assert A.norm().contains(a.norm())
            assert A.trace().contains(2*a.a)


@pytest.mark.parametrize("p", [2, 3, 5])
@pytest.mark.parametrize("n", [-3, 0, 1, 7, 19])
def test_absolute_precision_nonintegral_sources_and_reevaluation(p, n):
    place = GaussianLocalPlace(Qp(p, 3), root_sign=-1)
    source = GaussianBridge(F(1, p**6), F(-2, p**5))
    image = source.local(place, absprec=n)
    assert image.precision == n
    fine = source.local(place, absprec=25)
    assert fine.precision == 25 and fine.is_congruent(image)
    rational = RationalBridge(F(1, p**3)).local(place, absprec=n)
    assert rational.precision == n and rational.contains(F(1, p**3))
    assert source.local(place).precision == place.field.precision


def test_split_cancellation_requires_reevaluation_of_exact_source():
    place = GaussianLocalPlace(Qp(5, 3), root_sign=1)
    source = GaussianBridge(-7, 1)
    low = source.local(place, absprec=2)
    assert low.is_zero()
    with pytest.raises(PrecisionError):
        low.inverse()
    high = source.local(place, absprec=8)
    assert high.valuation() == 2
    assert (high*high.inverse()).is_congruent(1)
    assert source.local(place.conjugate(), absprec=8).valuation() == 0


@pytest.mark.parametrize("target", [
    Qp(5, 12), UnramifiedQuadraticField(Qp(3, 12), [1, 0, 1]),
    EisensteinQuadraticField(Qp(3, 12), [-3, 0, 1]),
    GaussianLocalPlace(Qp(2, 12), root_sign=-1),
])
def test_general_rational_quaternion_base_change(target):
    H = QuaternionAlgebra(a=2, b=-3)
    q, r = H(F(1, 3), 2, -1, 4), H(2, -3, 1, 1)
    Q, R = QuaternionBridge(q).local(target), QuaternionBridge(r).local(target)
    assert (Q.parent.a, Q.parent.b) == (2, -3)
    assert (Q*R).is_congruent(QuaternionBridge(q*r).local(target))
    assert Q.conjugate().is_congruent(QuaternionBridge(q.conjugate()).local(target))
    assert Q.norm().is_congruent(RationalBridge(q.norm()).local(target))
    assert Q.inverse().is_congruent(QuaternionBridge(q.inverse()).local(target))


@pytest.mark.parametrize("p", [2, 3, 5, 7, 13])
@pytest.mark.parametrize("sign", [-1, 1])
def test_biquaternion_operations_and_clifford_commuting_square(p, sign):
    place = GaussianLocalPlace(Qp(p, 14), root_sign=sign)
    q = Biquaternions(Gaussian(1, 2), Gaussian(3, -1), F(1, p), 1)
    r = Biquaternions(2, Gaussian(-1, 3), 4, Gaussian(1, -2))
    Q, R = QuaternionBridge(q).local(place), QuaternionBridge(r).local(place)
    assert (Q+R).is_congruent(QuaternionBridge(q+r).local(place))
    assert (Q*R).is_congruent(QuaternionBridge(q*r).local(place))
    assert Q.norm().is_congruent(place(q.norm()))
    assert Q.conjugate().is_congruent(QuaternionBridge(q.conjugate()).local(place))
    assert Q.inverse().is_congruent(QuaternionBridge(q.inverse()).local(place))
    C = CliffordAlgebra((1, 1, 1), place.field)
    mapped = q.to_clifford().map_coefficients(C, place)
    assert mapped.is_congruent(Q.to_clifford())
    assert (mapped*R.to_clifford()).is_congruent((Q*R).to_clifford())
    assert mapped.to_quaternion().is_congruent(Q)
    scalar_conjugate = Biquaternions(*(c.conjugate() for c in q.coefficients))
    assert QuaternionBridge(scalar_conjugate).local(place).is_congruent(
        QuaternionBridge(q).local(place.conjugate()))


@pytest.mark.parametrize("p", [2, 3, 5])
def test_independent_sympy_matrix_product_and_determinant(p):
    def expr(c):
        return sp.Rational(c.a.numerator, c.a.denominator)+sp.I*sp.Rational(c.b.numerator, c.b.denominator)

    def gaussian(z):
        re, im = sp.expand(z).as_real_imag()
        return Gaussian(F(re), F(im))

    def matrix(q):
        w, x, y, z = map(expr, q.coefficients)
        return sp.Matrix([[w+sp.I*x, y+sp.I*z], [-y+sp.I*z, w-sp.I*x]])

    q = Biquaternions(Gaussian(F(1, p), 2), Gaussian(1, -1), 3, Gaussian(-2, 1))
    r = Biquaternions(2, Gaussian(-1, 3), Gaussian(1, 2), -1)
    product = sp.expand(matrix(q)*matrix(r))
    w = (product[0, 0]+product[1, 1])/2
    x = (product[0, 0]-product[1, 1])/(2*sp.I)
    y = (product[0, 1]-product[1, 0])/2
    z = (product[0, 1]+product[1, 0])/(2*sp.I)
    place = GaussianLocalPlace(Qp(p, 16), root_sign=-1)
    computed = QuaternionBridge(q).local(place)*QuaternionBridge(r).local(place)
    for c, reference in zip(computed.coefficients, map(gaussian, (w, x, y, z))):
        assert c.is_congruent(place(reference))
    assert QuaternionBridge(q).local(place).norm().is_congruent(place(gaussian(matrix(q).det())))


@pytest.mark.parametrize("p", [2, 3, 5])
def test_sparse_clifford_exact_missing_blades_keep_requested_precision(p):
    place = GaussianLocalPlace(Qp(p, 3), root_sign=1)
    C, n = CliffordAlgebra((1, 1, 1), place.field), 17
    source = Hamilton(0, 2, 0, 0)
    mapped = source.to_clifford().map_coefficients(C, lambda c: place(c, absprec=n))
    assert mapped.to_quaternion().is_congruent(QuaternionBridge(source).local(place, absprec=n), n)
    assert all(c.precision == n for c in mapped.to_quaternion().coefficients)
    zero = C({6: place.field.O(n)})
    assert C({}).is_congruent(zero, n) and zero.is_congruent(C({}), n)
    assert (mapped-mapped).is_congruent(C({}), n)
    with pytest.raises(PrecisionError):
        C(0).is_congruent(zero, n)
    with pytest.raises(PrecisionError):
        C({}).is_congruent(C({6: place.field.O(n-1)}), n)


@pytest.mark.parametrize("p", [2, 3, 5])
def test_zero_reduced_norm_does_not_become_invertible(p):
    q = Biquaternions(1, Gaussian(0, 1), 0, 0)
    assert q.norm() == Gaussian(0)
    with pytest.raises(ZeroDivisionError):
        q.inverse()
    image = QuaternionBridge(q).local(GaussianLocalPlace(Qp(p), root_sign=1))
    assert image.norm().is_zero()
    with pytest.raises(PrecisionError):
        image.inverse()


def test_invalid_place_source_target_and_precision_are_rejected():
    with pytest.raises(TypeError):
        GaussianLocalPlace(Qp(3))
    with pytest.raises(TypeError):
        GaussianLocalPlace(3, root_sign=1)
    for sign in (True, 1.0, 0, 2, "1"):
        with pytest.raises(ValueError):
            GaussianLocalPlace(Qp(3), root_sign=sign)
    place = GaussianLocalPlace(Qp(3), root_sign=1)
    for source in (1.0, 1j, Qp(3)(1), QuadraticElement(Gaussian, Qp(3)(1), 0),
                   QuadraticElement(Gaussian, QuadraticElement(Gaussian, Qp(3)(1), 0), 0)):
        with pytest.raises(TypeError):
            GaussianBridge(source)
    with pytest.raises(TypeError):
        GaussianBridge(Gaussian(1), 2)
    for q in (Quaternion(Hamilton, (1, 2, 3)), Quaternion(Hamilton, (1.0, 2, 3, 4)),
              Quaternion(Biquaternions, (QuadraticElement(Gaussian, Qp(3)(1), 0), 0, 0, 0))):
        with pytest.raises(TypeError):
            QuaternionBridge(q)
    for target in (Qp(3), place.field):
        with pytest.raises(TypeError):
            GaussianBridge(1, 2).local(target)
        with pytest.raises(TypeError):
            QuaternionBridge(Biquaternions(1)).local(target)
    for target in (object(), Gaussian, QuadraticAlgebra(Qp(3))):
        with pytest.raises(TypeError):
            RationalBridge(1).local(target)
    for n in (True, 2.5, "3"):
        with pytest.raises(TypeError):
            place(1, absprec=n)
    with pytest.raises(DomainError):
        Hamilton(1).to_clifford().map_coefficients(CliffordAlgebra((1, 1, -1)))
    with pytest.raises(TypeError):
        Hamilton(1).to_clifford().map_coefficients(Qp(3))
