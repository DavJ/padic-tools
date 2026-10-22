from fractions import Fraction as F
from math import comb, factorial

import pytest

from padic_tools import (
    AtomicMeasure, DomainError, Polynomial, Qp, RationalAmiceMeasure,
    exp, valuation, zeta_special,
)


@pytest.mark.parametrize("p", [2,3,5,7])
def test_integer_atomic_conversion_preserves_transforms_and_moments(p):
    K = Qp(p,14)
    atom = AtomicMeasure(K,[(-3,F(1,p)),(0,7),(2,-3),(4,F(2,3))])
    measure = RationalAmiceMeasure.from_atomic(atom)
    assert measure.amice_series(14) == atom.amice_series(14)
    assert measure.laplace_series(14) == atom.laplace_series(14)
    for n in range(11):
        assert measure.moment(n) == atom.moment(n)
    for parameter in [p,p**3]:
        assert measure.amice(parameter).is_congruent(atom.amice(parameter))
    assert measure.laplace(p*p).is_congruent(atom.laplace(p*p))
    # An explicit high precision argument is not capped at the field default.
    value = measure.amice(K(p,absprec=30))
    assert value.precision >= 30+measure.valuation_bound


@pytest.mark.parametrize("p", [2,3,5,7])
def test_normalization_cancellation_and_proved_coefficient_bound(p):
    K = Qp(p)
    c = 3 if p == 2 else 2
    measure = RationalAmiceMeasure(K,[F(1,p),3],[c,1,2])
    assert measure.denominator[0] == 1
    assert all(valuation(a,p) >= measure.valuation_bound
               for a in measure.amice_series(60).coefficients)
    assert RationalAmiceMeasure(K,[p,1],[p,1]).moment(0) == 1
    removable = RationalAmiceMeasure(K,[0,2],[0,4])
    assert removable.moment(0) == F(1,2)
    assert removable.moment(3) == 0
    assert RationalAmiceMeasure(K,[1],[p]).valuation_bound == -1
    with pytest.raises(DomainError): RationalAmiceMeasure(K,[1],[p,1])
    with pytest.raises(DomainError): RationalAmiceMeasure(K,[1],[0,1])
    with pytest.raises(ZeroDivisionError): RationalAmiceMeasure(K,[1],[])
    zero = RationalAmiceMeasure(K,[],[p,1])
    assert zero.denominator == Polynomial([1])
    assert zero.moment(20) == 0
    assert zero.amice(p).contains(0)


@pytest.mark.parametrize("p", [3,5,7])
def test_rational_measure_exact_moments_match_zeta_special_values(p):
    # This is the complete rational transform, not a finite point measure.
    measure = RationalAmiceMeasure(Qp(p),[1],[2,1])
    assert measure.amice_series(30).coefficients == tuple(
        F((-1)**n,2**(n+1)) for n in range(30))
    for n in range(25):
        assert measure.moment(n) == (-1)**n*(1-2**(n+1))*zeta_special(-n)
    series = measure.laplace_series(16)
    assert series.coefficients == tuple((-1)**n*measure.moment(n)/factorial(n)
                                       for n in range(16))


@pytest.mark.parametrize("p,level", [(3,4),(5,3),(7,2)])
def test_rational_measure_matches_independent_finite_residue_sums(p,level):
    K = Qp(p,14)
    measure = RationalAmiceMeasure(K,1,[2,1])
    count = p**level
    # Cylinder values (-1)^a/2 are compatible when p is odd.
    for a in range(4):
        assert sum(F((-1)**(a+j*count),2) for j in range(p)) == F((-1)**a,2)
    for n in range(9):
        finite = sum((F((-1)**a,2)*a**n for a in range(count)),F(0))
        assert valuation(finite-measure.moment(n),p) >= level
    t = F(p)
    finite = sum((F((-1)**a,2)*(1+t)**a for a in range(count)),F(0))
    assert measure.amice(K(t)).is_congruent(finite,precision=level+1)


@pytest.mark.parametrize("p", [2,3,5])
def test_exact_polynomial_pairing_and_conditional_mahler_tail_enclosure(p):
    K = Qp(p,18)
    c = 3 if p == 2 else 2
    measure = RationalAmiceMeasure(K,[F(1,p)],[c,1])
    polynomial = Polynomial([F(1,7),-2,3,0,4])
    values = [polynomial(F(i)) for i in range(polynomial.degree+1)]
    differences = []
    while values:
        differences.append(values[0])
        values = [b-a for a,b in zip(values,values[1:])]
    expected = sum((v*measure.binomial_moment(n) for n,v in enumerate(differences)),F(0))
    assert measure.integrate_polynomial(polynomial) == expected
    # f(x)=(1+p)^x has exact Mahler coefficients p^n, hence tail >=N.
    for order in [1,4,10]:
        value = measure.integrate_mahler([p**n for n in range(order)],tail_valuation=order)
        assert value.precision == order+measure.valuation_bound
        assert value.is_congruent(measure.amice(p),precision=value.precision)
    assert measure.integrate_mahler([],tail_valuation=5).precision == 4
    with pytest.raises(TypeError): measure.integrate_mahler([K(1)],tail_valuation=3)
    with pytest.raises(TypeError): measure.integrate_mahler([1],tail_valuation=True)


@pytest.mark.parametrize("p", [2,3,5,7])
def test_measure_addition_scaling_and_additive_convolution(p):
    K = Qp(p,16)
    c = 3 if p == 2 else 2
    mu = RationalAmiceMeasure(K,[1],[c,1])
    nu = RationalAmiceMeasure(K,[1,2],[1,1])
    convolution = mu.convolve(nu)
    for n in range(9):
        assert convolution.moment(n) == sum(comb(n,k)*mu.moment(k)*nu.moment(n-k)
                                            for k in range(n+1))
        assert (mu+nu).moment(n) == mu.moment(n)+nu.moment(n)
        assert mu.scale(F(1,p)).moment(n) == mu.moment(n)/p
    t = K(p*p)
    assert convolution.laplace(t).is_congruent(mu.laplace(t)*nu.laplace(t))
    assert (mu+nu).amice(p).is_congruent(mu.amice(p)+nu.amice(p))


def test_domains_invalid_prefixes_and_unsupported_conversions():
    K = Qp(5)
    mu = RationalAmiceMeasure(K,1,[2,1])
    with pytest.raises(TypeError): RationalAmiceMeasure(5,1,[2,1])
    with pytest.raises(TypeError): RationalAmiceMeasure(K,1.0,[2,1])
    with pytest.raises(DomainError): mu.amice(1)
    with pytest.raises(DomainError): mu.laplace(1)
    with pytest.raises(TypeError): mu.amice(Qp(3)(3))
    with pytest.raises(TypeError): mu.convolve(RationalAmiceMeasure(Qp(3),1,[2,1]))
    with pytest.raises(TypeError): mu.integrate_polynomial([1,2])
    with pytest.raises(TypeError): RationalAmiceMeasure.from_atomic(mu)
    with pytest.raises(DomainError): RationalAmiceMeasure.from_atomic(AtomicMeasure(K,[(F(1,2),1)]))
    mu.amice_series(1)  # Cache entries must not bypass subsequent validation.
    for invalid in [True,-1,F(1,2)]:
        with pytest.raises(ValueError): mu.binomial_moment(invalid)
        with pytest.raises(ValueError): mu.moment(invalid)
        with pytest.raises(ValueError): mu.amice_series(invalid)
    with pytest.raises(ValueError): mu.amice_series(0)
    empty = RationalAmiceMeasure.from_atomic(AtomicMeasure(K))
    assert empty.moment(0) == 0
