from fractions import Fraction as F
from math import comb

import pytest

from padic_tools import AtomicMeasure, DomainError, Qp, PowerSeries


@pytest.mark.parametrize("p", [2, 3, 5, 7, 11])
def test_mellin_interpolates_exact_unit_moments_in_every_branch(p):
    K = Qp(p, 14)
    denominator = 3 if p == 2 else 2
    mu = AtomicMeasure(K, [(1, F(1,p)), (-1, -2),
                          (1+p, 3), (F(1,denominator), F(2,3))])
    period = 2 if p == 2 else p-1
    for n in range(2*period+2):
        value = mu.mellin(n, branch=n % period)
        assert value.precision >= 12
        assert value.contains(mu.moment(n))
    for n in [-1, -3]:
        value = mu.mellin(n, branch=n % period)
        assert value.contains(sum(w*a**n for a,w in mu.atoms))


@pytest.mark.parametrize("p", [2, 3, 5, 7])
def test_mellin_characters_are_multiplicative_at_noninteger_exponents(p):
    K = Qp(p, 15)
    s = F(1,3 if p == 2 else 2)
    a,b = 1+p, 1+2*p
    for branch in [0,1,-1]:
        product = (AtomicMeasure(K,[(a,1)]).mellin(s,branch=branch) *
                   AtomicMeasure(K,[(b,1)]).mellin(s,branch=branch))
        joint = AtomicMeasure(K,[(a*b,1)]).mellin(s,branch=branch)
        assert product.precision >= 14 and joint.precision >= 14
        assert product.is_congruent(joint)


def test_explicit_support_restriction_cancellation_and_validation():
    K = Qp(5)
    mu = AtomicMeasure(K,[(0,7),(5,2),(2,3),(2,-3),(1,4)])
    assert mu.atoms == ((F(0),F(7)),(F(1),F(4)),(F(5),F(2)))
    assert mu.restrict_units().atoms == ((F(1),F(4)),)
    with pytest.raises(DomainError): mu.mellin(0)
    assert mu.restrict_units().mellin(0).contains(4)
    with pytest.raises(DomainError): AtomicMeasure(K,[(F(1,5),1)])
    with pytest.raises(TypeError): AtomicMeasure(K,[(1.0,1)])
    with pytest.raises(TypeError): AtomicMeasure(K,[(1,K(2))])
    with pytest.raises(TypeError): mu.restrict_units().mellin(0,branch=True)
    with pytest.raises(DomainError): mu.restrict_units().mellin(F(1,5))
    with pytest.raises(TypeError): mu.laplace(Qp(3)(3))
    with pytest.raises(TypeError): mu.amice(Qp(3)(3))
    with pytest.raises(TypeError): mu.restrict_units().mellin(Qp(3)(3))
    for n in [-1, True, F(1,2)]:
        with pytest.raises(ValueError): mu.moment(n)


@pytest.mark.parametrize("p", [2, 3, 5, 7])
def test_laplace_additive_convolution_and_amice_relation(p):
    K = Qp(p, 16)
    atoms = [(F(-1),F(1,p)),(F(2),F(-3)),(F(0),F(7))]
    other = [(F(1),F(2)),(F(3),F(-1))]
    mu,nu = AtomicMeasure(K,atoms),AtomicMeasure(K,other)
    convolution = AtomicMeasure(K,[(a+b,w*v) for a,w in atoms for b,v in other])
    t = K(p*p)
    value = convolution.laplace(t)
    assert value.precision >= 14
    assert value.is_congruent(mu.laplace(t)*nu.laplace(t))
    # Exact coefficient identity A_mu(exp(-t)-1)=L_mu(t).
    z = PowerSeries([0,1],order=13)
    assert mu.amice_series(13).compose((-z).exp()-1) == mu.laplace_series(13)
    from padic_tools import exp
    assert mu.amice(exp(-t)-1).is_congruent(mu.laplace(t))


def test_exact_laplace_taylor_moments_and_amice_binomial_moments():
    mu = AtomicMeasure(Qp(5),[(0,1),(2,3),(4,-2)])
    assert mu.laplace_series(4).coefficients == (2,2,-10,F(52,3))
    assert mu.amice_series(6).coefficients == tuple(
        (1 if n == 0 else 0)+3*comb(2,n)-2*comb(4,n)
        for n in range(6))
    for method in [mu.amice_series,mu.laplace_series]:
        for n in [0,-1,True]:
            with pytest.raises(ValueError): method(n)


@pytest.mark.parametrize("p,point", [(2,F(1,3)),(3,F(1,2)),(5,F(-2,3))])
def test_amice_binomial_evaluation_on_entire_open_unit_disc(p,point):
    K = Qp(p, 18)
    mu = AtomicMeasure(K,[(point,F(1,p))])
    t = K(p,absprec=17)
    value = p*mu.amice(t)
    assert value.precision >= 17
    # Rational exponent is verified by an integer-power identity.
    assert (value**point.denominator).is_congruent((1+t)**point.numerator)
    with pytest.raises(DomainError): mu.amice(K(1))


def test_transform_domains_use_support_and_dyadic_convergence():
    for p in [2,3,5]:
        K = Qp(p)
        mu = AtomicMeasure(K,[(1,1)])
        with pytest.raises(DomainError): mu.laplace(K(1))
        if p == 2:
            with pytest.raises(DomainError): mu.laplace(K(2))
        assert mu.laplace(K(p*p)).precision >= K.precision
        assert AtomicMeasure(K,[(p*p,1)]).laplace(K(1)).precision >= K.precision
        # Coarse balls cannot certify the local exponential domain.
        with pytest.raises(DomainError): mu.laplace(K.O(0))
        empty = AtomicMeasure(K)
        assert empty.moment(0) == 0
        assert empty.laplace(1).contains(0)
        assert empty.mellin(0).contains(0)
        assert empty.amice(K(p)).contains(0)
        assert AtomicMeasure(K,[(0,3)]).laplace(1).contains(3)
