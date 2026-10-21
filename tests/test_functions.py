from fractions import Fraction as F
from math import factorial
import pytest

from padic_tools import (
    Qp, DomainError, PrecisionError, exp, log, log1p, sin, cos, sqrt, power,
    teichmuller, theta3, theta4, zeta_special, zeta_euler_value, zeta_p, bernoulli,
)


@pytest.mark.parametrize("p", [2,3,5,7,11])
def test_analytic_identities_and_stability(p):
    K=Qp(p,18); step=4 if p==2 else p
    x,y=K(step),K(2*step)
    assert log(exp(x)).is_congruent(x)
    assert exp(log(1+x)).is_congruent(1+x)
    assert exp(x+y).is_congruent(exp(x)*exp(y))
    assert log((1+x)*(1+y)).is_congruent(log(1+x)+log(1+y))
    assert (sin(x)**2+cos(x)**2).is_congruent(1)
    assert power(1+x,3).is_congruent((1+x)**3)
    high=Qp(p,28)
    assert exp(high(step)).truncate(18).is_congruent(exp(x))
    assert log1p(high(step)).truncate(18).is_congruent(log1p(x))
    # An independent exact rational partial sum with a much longer tail.
    oracle=sum((F(step**k,factorial(k)) for k in range(100)),F(0))
    assert exp(x).contains(oracle)


def test_domain_boundaries_and_log_branch():
    with pytest.raises(DomainError): exp(Qp(2)(2))
    with pytest.raises(DomainError): exp(Qp(5)(1))
    with pytest.raises(DomainError): log1p(Qp(5)(1))
    with pytest.raises(PrecisionError): log(Qp(3).O(8))
    assert log(Qp(5)(5)).is_zero()
    assert log(Qp(5)(5),p_branch=7).is_congruent(7)
    assert log1p(Qp(2)(-2)).is_zero()  # log(-1)=0, including denominator spikes.
    with pytest.raises(PrecisionError): exp(Qp(5)(5),precision=100)
    assert exp(Qp(3).O(10)).is_congruent(1)


@pytest.mark.parametrize("p", [2,3,5,7,11,13,17,101])
def test_square_root_for_valuations_and_precision(p):
    K=Qp(p,14)
    for q in (F(1),F(9),F(25),F(1,p*p),F(4*p*p)):
        x=K(q); y=sqrt(x)
        assert (y*y).is_congruent(x)
        expected=x.precision-x.v//2-(1 if p==2 else 0)
        assert y.precision==expected
    if p!=2:
        for a in range(1,min(p,20)):
            if pow(a,(p-1)//2,p)==p-1:
                with pytest.raises(DomainError): sqrt(K(a))
    with pytest.raises(DomainError): sqrt(K(p))
    with pytest.raises(PrecisionError): sqrt(K.O(12))


def test_dyadic_square_root_boundaries():
    K=Qp(2)
    with pytest.raises(DomainError): sqrt(K(3))
    with pytest.raises(PrecisionError): sqrt(K(1,absprec=2))
    assert sqrt(K(9)).residue(2)==1


@pytest.mark.parametrize("p", [2,3,5,7,11])
def test_teichmuller_and_unit_decomposition(p):
    K=Qp(p,15)
    for a in range(1,12):
        if a%p:
            w=teichmuller(K(a))
            assert (w**(2 if p==2 else p-1)).is_congruent(1)
            assert log(w).is_zero()
            assert exp(log(K(a))).is_congruent(K(a)/w)


@pytest.mark.parametrize("p", [2,3,5,7])
def test_theta_against_finite_exact_sum_and_high_precision(p):
    K=Qp(p,16)
    for z in (F(1),F(-1),F(p),F(1,p*p)):
        q=K(p)
        value=theta3(q,z)
        oracle=1+sum((F(p)**(n*n)*(z**n+z**(-n)) for n in range(1,20)),F(0))
        assert value.contains(oracle)
        assert value.is_congruent(theta3(Qp(p,24)(p),z))
    assert theta4(K(p)).is_congruent(theta3(K(p),-1))
    assert theta3(K.O(16)).is_congruent(1)
    with pytest.raises(DomainError): theta3(K(1))
    with pytest.raises(DomainError): theta3(K(p),0)


@pytest.mark.parametrize("p", [2,3,5,7,11])
@pytest.mark.parametrize("n", [2,4,6,8,10,12,16,24])
def test_zeta_branches_against_exact_interpolation(p,n):
    K=Qp(p,12)
    r=0 if p==2 else n%(p-1)
    assert zeta_p(1-n,K,branch=r).is_congruent(zeta_euler_value(1-n,p))


@pytest.mark.parametrize("p", [2,3,5,7])
def test_zeta_at_positive_rational_and_finite_input(p):
    K=Qp(p,12)
    for s in (F(2),F(3),F(1,1+p)):
        low=zeta_p(s,K)
        high=zeta_p(s,Qp(p,20))
        assert low.precision==12 and low.is_congruent(high)
        uncertain=zeta_p(K(s))
        assert uncertain.is_congruent(high)
    with pytest.raises(DomainError): zeta_p(1,K)
    with pytest.raises(PrecisionError): zeta_p(K(1))
    with pytest.raises(DomainError): zeta_p(F(1,p),K)
    with pytest.raises(DomainError): zeta_p(-1,K,branch=1)


def test_zeta_wrong_branch_is_not_silently_a_special_value():
    K=Qp(5,12)
    assert not zeta_p(-1,K,branch=0).is_congruent(zeta_euler_value(-1,5))
    assert zeta_special(0)==F(-1,2)
    assert zeta_special(-1)==F(-1,12)
    assert zeta_special(-3)==F(1,120)
    assert zeta_special(-2)==0
    assert bernoulli(1)==F(-1,2) and bernoulli(12)==F(-691,2730)
