from fractions import Fraction as F
import pytest
from padic_tools import (
    Polynomial as P, PowerSeries as S, Qp, hensel_lift, DomainError,
    exp_series, log1p_series, theta3_series, exp, zeta_p,
)


def test_polynomial_euclidean_arithmetic_and_calculus():
    x=P([0,1]); f=(x-1)**3*(x+2)
    q,r=divmod(f,(x-1)**2)
    assert r==P() and q==(x-1)*(x+2)
    assert f.gcd(f.derivative())==(x-1)**2
    assert f.derivative().integral(f[0])==f
    assert (x*x-2)(F(3))==7
    with pytest.raises(ZeroDivisionError): divmod(f,P())


@pytest.mark.parametrize("p,root", [(3,1),(5,1),(7,1),(11,1)])
def test_hensel_lifts_exact_polynomial_roots(p,root):
    f=P([-1-p,0,1]); K=Qp(p,12)
    r=hensel_lift(f,root,K)
    assert f(r).is_zero()
    assert r.residue(1)==root
    with pytest.raises(DomainError): hensel_lift(P([0,0,1]),0,K)
    with pytest.raises(DomainError): hensel_lift(P([F(1,p),1]),0,K)


def test_formal_series_and_order_semantics():
    x=S([0,1],12)
    assert x.exp()==exp_series(12)
    assert (1+x).log()==log1p_series(12)
    assert x.exp().log()==x
    assert ((1+x)*(1+x).inverse())==S([1],12)
    assert x.exp().compose(x*x)==(x*x).exp()
    assert x.derivative().order==11
    assert x.integral().order==13
    assert (x+S([0,1],7)).order==7
    assert theta3_series(12).coefficients==(1,2,0,0,2,0,0,0,0,2,0,0)
    assert S([1],1).log()==S([],1)
    with pytest.raises(DomainError): (1+x).exp()
    with pytest.raises(DomainError): x.log()
    with pytest.raises(DomainError): x.compose(1+x)
    with pytest.raises(ZeroDivisionError): x.inverse()


def test_sympy_conversion_and_function_evaluation():
    sp=pytest.importorskip("sympy")
    from padic_tools.sympy_bridge import PadicExp,PadicLog,PadicTheta3,PadicZeta,evaluate
    t=sp.Symbol("t")
    assert P([F(1,3),0,2]).to_sympy("t")==sp.Rational(1,3)+2*t*t
    assert exp_series(8).to_sympy("t").getO()==sp.Order(t**8,(t,0))
    assert sp.diff(PadicExp(t),t)==PadicExp(t)
    assert sp.diff(PadicLog(t),t)==1/t
    K=Qp(5,15)
    assert evaluate(PadicExp(t),K,{t:5}).is_congruent(exp(K(5)))
    assert evaluate(PadicLog(PadicExp(t)),K,{t:5}).is_congruent(5)
    assert evaluate(PadicTheta3(t,1),K,{t:5}).precision==15
    assert evaluate(PadicZeta(-3,5,0),K).is_congruent(zeta_p(-3,K))
    assert evaluate((t+1)**2,K,{t:F(1,3)}).contains(F(16,9))
    with pytest.raises(DomainError): evaluate(sp.zeta(t),K,{t:2})
    with pytest.raises(DomainError): evaluate(sp.sqrt(t),K,{t:1})
    with pytest.raises(DomainError): evaluate(PadicZeta(-3,3,0),K)
    with pytest.raises(ValueError): evaluate(t,K)
