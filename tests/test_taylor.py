from fractions import Fraction as F
from math import factorial

import pytest

from padic_tools import (
    DomainError, Polynomial, PowerSeries, Qp, TaylorExpansion,
    cos, exp, exp_series, log, log1p, log1p_series, sin, taylor,
)


def test_exact_polynomial_translation_and_taylor_derivative_identity():
    x = Polynomial([0,1])
    f = (x*x-3*x+F(2,7))**4
    center = F(2,3)
    shifted = f.shift(center)
    for h in [F(0),F(1),F(-2,5)]:
        assert shifted(h) == f(center+h)
    jet = f.taylor(center,order=12)
    derivative = f
    for n,c in enumerate(jet.coefficients):
        assert c == derivative(center)/factorial(n)
        derivative = derivative.derivative()
    assert Polynomial().shift(center) == Polynomial()
    assert Polynomial().taylor(center,order=2) == PowerSeries([],2)


def test_formal_standard_taylor_jets_and_supported_functions():
    assert taylor(exp,order=13) == exp_series(13)
    assert taylor(log1p,order=13) == log1p_series(13)
    assert taylor(log,1,order=13) == log1p_series(13)
    s,c = taylor(sin,order=15),taylor(cos,order=15)
    assert s.derivative() == PowerSeries(c.coefficients,14)
    assert c.derivative() == PowerSeries((-s).coefficients,14)
    assert s*s+c*c == PowerSeries([1],15)
    with pytest.raises(TypeError): taylor(exp,1)
    with pytest.raises(TypeError): taylor(lambda x:x*x,Qp(5)(1))
    with pytest.raises(TypeError): TaylorExpansion()
    for order in [0,-1,True,F(1,2)]:
        with pytest.raises(ValueError): taylor(exp,order=order)


@pytest.mark.parametrize("p", [2,3,5,7])
@pytest.mark.parametrize("function", [exp,sin,cos,log1p])
def test_local_analytic_taylor_encloses_independent_high_precision_values(p,function):
    K,oracle = Qp(p,18),Qp(p,35)
    center = p*p
    point = center+p**3
    jet = taylor(function,K(center),order=30)
    value = jet(K(point))
    assert value.precision >= 18
    assert value.is_congruent(function(oracle(point)),precision=value.precision)
    # Every coefficient agrees with the exact derivative formula for exp.
    if function is exp:
        for n,c in enumerate(jet.coefficients):
            assert (c*factorial(n)).is_congruent(exp(oracle(center)))


@pytest.mark.parametrize("p", [2,3,5])
@pytest.mark.parametrize("scale_power", [-2,0,2])
def test_log_taylor_at_nonunit_centers_and_its_local_domain(p,scale_power):
    K = Qp(p,22)
    center = F(p)**scale_power * (1+p*p)
    point = center*(1+p*p)
    jet = taylor(log,K(center),order=20)
    value = jet(K(point))
    assert value.precision >= 18
    assert value.is_congruent(log(Qp(p,40)(point)),precision=value.precision)
    with pytest.raises(DomainError): jet(K(center*p))


@pytest.mark.parametrize("p", [2,3,5])
def test_taylor_remainder_includes_later_denominator_spikes(p):
    K = Qp(p,30)
    h = p*p if p == 2 else p
    for order in range(1,12):
        for function in [exp,log1p,sin,cos]:
            jet = taylor(function,K(0),order=order)
            value = jet(K(h))
            assert value.precision >= jet.tail_precision(K(h))
            assert value.contains(function(Qp(p,45)(h)).representative())


def test_polynomial_taylor_remainder_nonintegral_points_and_short_orders():
    K = Qp(5,20)
    f = Polynomial([F(1,25),-7,0,2,4])
    for order in [1,2,5,8]:
        jet = f.taylor(K(F(1,5)),order=order)
        for point in [F(2),F(-2,25),F(1,5)]:
            assert jet(K(point)).contains(f(point))
        if order >= 5:
            assert jet.tail_precision(K(2)) is None
    assert taylor(Polynomial(),K(2),order=3)(K(7)).contains(0)


def test_taylor_convergence_and_prime_checks_cover_entire_balls():
    K = Qp(5)
    for function in [exp,sin,cos,log1p]:
        with pytest.raises(DomainError): taylor(function,K(1))
    with pytest.raises(ArithmeticError): taylor(log,K.O(2))
    jet = taylor(exp,K(5),order=10)
    with pytest.raises(DomainError): jet(K(1))
    with pytest.raises(DomainError): jet(K.O(0))
    with pytest.raises(TypeError): jet(Qp(3)(3))
    # The 2-adic exponential requires a displacement divisible by 4.
    with pytest.raises(DomainError): taylor(exp,Qp(2)(4),order=10)(Qp(2)(6))
