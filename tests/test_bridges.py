from fractions import Fraction as F
import random
import pytest

from padic_tools import (
    Qp, DomainError, PrecisionError, Gaussian, Hamilton, Biquaternions,
    CliffordAlgebra, QuadraticAlgebra, QuaternionAlgebra, RationalBridge,
    GaussianBridge, QuaternionBridge, riemann_zeta,
)


def test_gaussian_exact_algebra_and_zero_divisor():
    i=Gaussian.generator(); z=Gaussian(F(1,3),F(2,5))
    assert i*i==Gaussian(-1)
    assert z*z.inverse()==Gaussian(1)
    assert z*z.conjugate()==Gaussian(z.norm())
    split=QuadraticAlgebra(d=1)
    assert split(1,1)*split(1,-1)==split(0)
    with pytest.raises(ZeroDivisionError): split(1,1).inverse()


@pytest.mark.parametrize("p", [2,3,5,7,13])
def test_diagonal_gaussian_embedding_preserves_rational_algebra(p):
    K=Qp(p,16); a=Gaussian(F(1,3),F(2,7)); b=Gaussian(F(-5,2),F(3,11))
    A=GaussianBridge(a).padic_algebra(K); B=GaussianBridge(b).padic_algebra(K)
    assert (A+B).is_congruent(GaussianBridge(a+b).padic_algebra(K))
    assert (A*B).is_congruent(GaussianBridge(a*b).padic_algebra(K))
    assert A.inverse().is_congruent(GaussianBridge(a.inverse()).padic_algebra(K))
    assert RationalBridge(F(1,3)).padic(K).contains(F(1,3))


@pytest.mark.parametrize("p", [5,13,17])
def test_split_places_are_explicit_and_distinct(p):
    K=Qp(p,15)
    bridge=GaussianBridge(0,1)
    plus,minus=bridge.split_padic(K),bridge.split_padic(K,root_sign=-1)
    assert (plus*plus).is_congruent(-1)
    assert plus.is_congruent(-minus) and not plus.is_congruent(minus)
    a,b=Gaussian(F(1,p),2),Gaussian(3,F(2,p))
    assert (GaussianBridge(a).split_padic(K)*GaussianBridge(b).split_padic(K)).is_congruent(
        GaussianBridge(a*b).split_padic(K))


def test_nonsplit_gaussian_and_padic_zero_divisors():
    for p in (2,3,7):
        with pytest.raises(DomainError): GaussianBridge(0,1).split_padic(Qp(p))
    K=Qp(5,12)
    root=GaussianBridge(0,1).split_padic(K)
    A=QuadraticAlgebra(K)
    divisor=A(root,1)
    # In the split algebra (r+i)(r-i)=0 because r^2=-1.
    assert divisor.norm().is_zero()
    with pytest.raises(PrecisionError): divisor.inverse()
    with pytest.raises(TypeError): GaussianBridge(A(1,2))


@pytest.mark.parametrize("parameters", [(-1,-1),(2,3),(-2,5)])
def test_general_quaternion_algebra(parameters):
    H=QuaternionAlgebra(a=parameters[0],b=parameters[1])
    i,j,k=H.generators()
    assert i*i==H(parameters[0]) and j*j==H(parameters[1])
    assert i*j==k and j*i==-k
    a,b=H(1,2,3,4),H(4,3,2,1)
    assert (a*b).norm()==a.norm()*b.norm()
    assert a*a.conjugate()==H(a.norm())
    assert a*a.inverse()==H(1)


def test_clifford_associativity_grades_and_involutions():
    C=CliffordAlgebra((1,-1,2))
    vectors=[C.e(i) for i in range(3)]
    for i,e in enumerate(vectors):
        assert e*e==C(C.metric[i])
        assert e^e==C(0)
        for j,f in enumerate(vectors):
            if i!=j: assert e*f==-f*e
    blades=[C({mask:1}) for mask in range(8)]
    for a in blades:
        for b in blades:
            for c in blades:
                assert (a*b)*c==a*(b*c)
    x=C({0:2,1:3,3:4,7:5}); y=C({2:-1,5:2})
    assert (x*y).reverse()==y.reverse()*x.reverse()
    assert x.involute().involute()==x
    assert x.conjugate().conjugate()==x
    assert x.grade(2)==C({3:4})


def test_exact_quaternion_clifford_isomorphism():
    rng=random.Random(124)
    for H in (Hamilton,Biquaternions):
        for _ in range(20):
            a=H(*(rng.randint(-5,5) for _ in range(4)))
            b=H(*(rng.randint(-5,5) for _ in range(4)))
            assert a.to_clifford().to_quaternion()==a
            assert (a*b).to_clifford()==a.to_clifford()*b.to_clifford()
            assert a.conjugate().to_clifford()==a.to_clifford().reverse()
    i,j,k=Hamilton.generators()
    assert i.to_clifford()*j.to_clifford()==k.to_clifford()
    C=CliffordAlgebra()
    assert C.pseudoscalar()*C.pseudoscalar()==C(-1)
    with pytest.raises(DomainError): C.e(0).to_quaternion()


@pytest.mark.parametrize("p", [2,3,5,7])
def test_biquaternion_padic_base_change_and_clifford_commute(p):
    K=Qp(p,14)
    a=Biquaternions(Gaussian(1,2),Gaussian(3,-4),Gaussian(F(2,7),5),1)
    b=Biquaternions(Gaussian(5,1),Gaussian(-2,3),1,Gaussian(0,2))
    A,B=QuaternionBridge(a).padic(K),QuaternionBridge(b).padic(K)
    assert (A*B).is_congruent(QuaternionBridge(a*b).padic(K))
    assert A.inverse().is_congruent(QuaternionBridge(a.inverse()).padic(K))
    assert (A*B).to_clifford().is_congruent(A.to_clifford()*B.to_clifford())
    assert A.to_clifford().to_quaternion().is_congruent(A)


def test_complex_images_and_classical_zeta():
    mp=pytest.importorskip("mpmath")
    before=mp.mp.dps
    z=GaussianBridge(F(1,3),F(2,7)).complex(70)
    assert abs(z.real*3-1)<mp.mpf("1e-65")
    assert mp.mp.dps==before
    assert abs(riemann_zeta(-3)-F(1,120))<mp.mpf("1e-48")
    for H in (Hamilton,Biquaternions):
        a=H(1,2,3,4)
        A=QuaternionBridge(a).complex(60)
        expected=QuaternionBridge(a*a).complex(60)
        assert all(abs(x-y)<mp.mpf("1e-55") for x,y in zip((A*A).coefficients,expected.coefficients))
    with pytest.raises(TypeError): riemann_zeta(Qp(5)(2))
    with pytest.raises(TypeError): riemann_zeta(0.5)
