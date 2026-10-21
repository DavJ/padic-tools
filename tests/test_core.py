from decimal import Decimal
from fractions import Fraction as F
import random

import pytest
from padic_tools import Qp, Padic, PrecisionError, DomainError, valuation
from padic_tools.core import is_prime


@pytest.mark.parametrize("p", [2,3,5,7,11,101])
def test_exact_rationals_and_residues(p):
    K=Qp(p,12)
    for n in range(-9,10):
        for d in (1,2,3,5,7,25):
            q=F(n,d); x=K(q)
            assert x.contains(q)
            if n:
                assert x.valuation()==valuation(q,p)
                assert x.relative_precision==12
                assert x.inverse().contains(1/q)
            if x.v>=0:
                assert valuation(F(x.residue())-q,p)>=x.precision


@pytest.mark.parametrize("p", [2,3,5,7])
def test_arithmetic_encloses_exact_rational_oracle(p):
    rng=random.Random(91283)
    K=Qp(p,9)
    for _ in range(100):
        a=F(rng.randint(-200,200),rng.randint(1,50))
        b=F(rng.randint(1,200),rng.randint(1,50))
        x,y=K(a),K(b)
        assert (x+y).contains(a+b)
        assert (x-y).contains(a-b)
        assert (x*y).contains(a*b)
        assert (x/y).contains(a/b)
        assert (x+b).contains(a+b)
        assert (x*b).contains(a*b)
        assert (x/b).contains(a/b)
        for n in (-2,0,1,3):
            if n>=0 or a:
                assert (x**n).contains(a**n)


@pytest.mark.parametrize("p", [2,3,5])
def test_ball_error_propagation_for_perturbed_centres(p):
    K=Qp(p,10)
    x=K(F(1,p),absprec=4); y=K(p,absprec=6)
    for k in (-2,-1,0,1,2):
        a=x.representative()+k*p**4
        b=y.representative()+(k+1)*p**6
        assert (x+y).contains(a+b)
        assert (x*y).contains(a*b)
        assert (x/y).contains(a/b)


def test_precision_formulas_and_scalar_identity():
    K=Qp(5,12); x=K(25,absprec=10); y=K(F(1,5),absprec=8)
    assert (x+y).precision==8
    assert (x*y).precision==9
    assert x.inverse().precision==6
    assert (x*25).precision==12
    assert (x/25).precision==8
    assert x+0==x and x*1==x and x/1==x
    assert (x-x).is_zero()
    assert (x-x).precision==10
    assert x.truncate(7).precision==7
    with pytest.raises(PrecisionError): x.truncate(11)
    with pytest.raises(PrecisionError): K(x,absprec=11)


def test_zero_balls_negative_precision_and_congruence():
    K=Qp(3)
    z=K.O(-2)
    assert z.contains(F(1,9)) and not z.contains(F(1,27))
    with pytest.raises(PrecisionError): z.valuation()
    with pytest.raises(PrecisionError): z.inverse()
    with pytest.raises(TypeError): bool(z)
    a,b=K(1,absprec=2),K(10,absprec=4)
    assert a!=b and a.is_congruent(b)
    assert len({a,b})==2
    with pytest.raises(PrecisionError): a.is_congruent(b,3)


def test_float_rejection_and_prime_validation():
    for q in (0,1,4,9,341,561,3215031751):
        assert not is_prime(q)
        with pytest.raises(ValueError): Qp(q)
    assert is_prime(2**61-1)
    with pytest.raises(ValueError): Qp(2**89-1)
    for x in (0.1,1+2j,True):
        with pytest.raises(TypeError): Qp(5)(x)
    assert Qp(5)(Decimal("0.125")).contains(F(1,8))
    assert Qp(5)("1/3").contains(F(1,3))
    with pytest.raises(TypeError): Qp(3)(1)+Qp(5)(1)
    with pytest.raises(DomainError): Qp(5)(F(1,5)).residue()


def test_no_accidental_public_uninitialized_constructor():
    with pytest.raises(TypeError): Padic(5,3)
