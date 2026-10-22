from fractions import Fraction as F
import random

import pytest

from padic_tools import (
    CyclotomicValue as CV, DomainError, PadicTestFunction as TF,
    PrecisionError, Qp, additive_character as chi, additive_phase,
    finite_fourier, haar_volume,
)


@pytest.mark.parametrize("p", [2,3,5,7])
def test_exact_character_is_additive_and_trivial_exactly_on_integral_arguments(p):
    K = Qp(p,12)
    points = [F(-2,p*p),F(1,p),F(5,p*p*(p+1)),F(1,p+1),F(0),F(7)]
    for x in points:
        phase = additive_phase(x,K)
        assert 0 <= phase < 1
        assert (chi(x,K) == 1) == (x.denominator % p != 0)
        assert chi(K(x)) == chi(x,p)
        assert chi(x,p)*chi(-x,p) == 1
        for y in points:
            assert chi(x+y,p) == chi(x,p)*chi(y,p)
    assert additive_phase(F(1,p*(p+1)),p) == F(1,p)
    assert additive_phase(F(-1,p),p) == F(p-1,p)
    assert chi(K(F(1,p),absprec=0)) == chi(F(1,p),p)
    with pytest.raises(PrecisionError): chi(K.O(-1))
    with pytest.raises(TypeError): chi(K(1),Qp(11))
    with pytest.raises(TypeError): chi(F(1,p))
    with pytest.raises(TypeError): additive_phase(0.2,p)


@pytest.mark.parametrize("p,level", [(2,1),(2,4),(3,1),(3,3),(5,2),(7,1)])
def test_cyclotomic_reduction_conjugation_and_level_descent(p,level):
    root = CV.root(p,level)
    n = p**level
    assert root**n == 1
    assert root.conjugate()*root == 1
    assert sum((CV.root(p,level,k) for k in range(n)),CV(p)) == 0
    if level > 1:
        assert root**p == CV.root(p,level-1)
        assert hash(root**p) == hash(CV.root(p,level-1))
    assert CV(p,[F(2,3)]) == F(2,3)
    assert hash(CV(p,[F(2,3)])) == hash(F(2,3))
    assert bool(root) and not CV(p)
    assert (3*root/2+root/2) == 2*root
    assert root.conjugate().conjugate() == root
    assert CV(p,[1,-1],level=0).as_fraction() == 0
    if n > 2:
        with pytest.raises(DomainError): root.as_fraction()
        with pytest.raises(DomainError): root/root


@pytest.mark.parametrize("p,level", [(2,3),(3,2),(5,1)])
def test_cyclotomic_arithmetic_against_independent_sympy_polynomial_quotient(p,level):
    sp = pytest.importorskip('sympy')
    rng = random.Random(430+p)
    x = sp.Symbol('x')
    n = p**level
    modulus = sp.cyclotomic_poly(n,x)
    def polynomial(value):
        step = p**(level-value.level)
        return sum(sp.Rational(c.numerator,c.denominator)*x**(i*step)
                   for i,c in enumerate(value.coefficients))
    for _ in range(8):
        ac = [F(rng.randint(-4,4),rng.randint(1,5)) for _ in range(n+2)]
        bc = [F(rng.randint(-4,4),rng.randint(1,5)) for _ in range(n)]
        a,b = CV(p,ac,level=level),CV(p,bc,level=level)
        original = sum(sp.Rational(c.numerator,c.denominator)*x**i for i,c in enumerate(ac))
        assert sp.rem(polynomial(a)-original,modulus,x) == 0
        assert sp.rem(polynomial(a+b)-polynomial(a)-polynomial(b),modulus,x) == 0
        assert sp.rem(polynomial(a*b)-polynomial(a)*polynomial(b),modulus,x) == 0


@pytest.mark.parametrize("p,level", [(2,0),(2,2),(2,4),(3,1),(3,2),(5,1),(7,1)])
def test_finite_transform_inverse_reflection_parseval_and_modulation(p,level):
    n = p**level
    values = tuple(F(j*j-2*j+3,j+1) for j in range(n))
    spectrum = finite_fourier(values,p)
    assert finite_fourier(spectrum,p,inverse=True) == values
    assert finite_fourier(spectrum,p) == tuple(F(values[-j % n],n) for j in range(n))
    assert sum(v*v.conjugate() for v in spectrum) == sum(v*v for v in values)/n
    mode = tuple(CV.root(p,level,3*j) for j in range(n))
    transformed = finite_fourier(mode,p)
    assert transformed == tuple(int(k == 3 % n) for k in range(n))
    assert spectrum[0] == sum(values)/n


@pytest.mark.parametrize("p,a,b", [(2,-2,1),(2,1,3),(3,-1,1),(3,0,2),(5,0,1)])
def test_qp_fourier_has_haar_scaling_dual_lattice_and_exact_inversion(p,a,b):
    K = Qp(p,10)
    n = p**(b-a)
    f = TF(K,[F(j-2,j+1) for j in range(n)],support_exponent=a,constancy_exponent=b)
    transformed = f.fourier()
    assert (transformed.support_exponent,transformed.constancy_exponent) == (-b,-a)
    assert transformed.fourier(inverse=True) == f
    assert transformed.fourier() == f.reflect()
    assert transformed.norm_squared() == f.norm_squared()
    assert transformed(0) == f.integral()
    assert transformed.integral() == f(0)
    step = F(p**a) if a >= 0 else F(1,p**(-a))
    for k in range(n):
        frequency = F(k,p**b) if b >= 0 else k*p**(-b)
        expected = haar_volume(p,b)*sum((f.values[j]*chi(-step*j*frequency,p)
                                         for j in range(n)),CV(p))
        assert transformed(frequency) == expected
    translated = f.translate(step)
    for k in range(n):
        frequency = F(k,p**b) if b >= 0 else k*p**(-b)
        assert translated.fourier().values[k] == transformed.values[k]*chi(-step*frequency,p)


@pytest.mark.parametrize("p", [2,3,5])
def test_haar_convolution_theorem_and_cross_parseval(p):
    K = Qp(p)
    f = TF(K,[F(j+1,p) for j in range(p)],support_exponent=-1,constancy_exponent=0)
    g = TF(K,[CV.root(p,1,j)+j for j in range(p)],support_exponent=-1,constancy_exponent=0)
    convolved = f.convolve(g)
    assert convolved.fourier() == f.fourier()*g.fourier()
    assert f.inner_product(g) == f.fourier().inner_product(g.fourier())
    assert convolved.integral() == f.integral()*g.integral()
    assert f.convolve(g) == g.convolve(f)
    assert (f+f)*2 == f*4
    assert CV.root(p,1)*f == f*CV.root(p,1)


@pytest.mark.parametrize("p,e", [(2,-2),(2,0),(2,3),(3,-1),(3,2),(5,1)])
def test_ball_volumes_and_analytic_fourier_indicator_formula(p,e):
    K = Qp(p,10)
    ball = TF.ball(K,0,e)
    expected = TF.ball(K,0,-e)*haar_volume(K,e)
    assert ball.fourier() == expected
    assert ball.integral() == haar_volume(p,e)
    assert ball.norm_squared() == haar_volume(p,e)
    assert ball.convolve(ball) == ball*haar_volume(p,e)
    unit_ball = TF.ball(K)
    assert unit_ball.fourier() == unit_ball


@pytest.mark.parametrize("p", [2,3,5])
def test_noncentral_ball_evaluation_and_its_fourier_phase(p):
    K = Qp(p,10)
    center,e = F(1,p),1
    f = TF.ball(K,center,e)
    assert f(center) == 1 and f(center+p**e) == 1
    assert f(0) == 0 and f(F(1,p*p)) == 0
    assert f(K(center,absprec=e)) == 1
    with pytest.raises(PrecisionError): f(K(center,absprec=e-1))
    assert TF.ball(K,K(center,absprec=e),e) == f
    with pytest.raises(PrecisionError): TF.ball(K,K.O(e-1),e)
    transformed = f.fourier()
    for k in range(len(transformed.values)):
        xi = F(k,p)
        assert transformed(xi) == haar_volume(p,e)*chi(-center*xi,p)
    assert transformed(F(1,p*p)) == 0
    assert f(K(F(1,p*p),absprec=-1)) == 0


def test_complex_and_symbolic_views_preserve_global_precision():
    mp = pytest.importorskip('mpmath')
    sp = pytest.importorskip('sympy')
    before = mp.mp.dps
    root = CV.root(2,2)
    assert root.to_sympy() == sp.I
    assert abs(root.complex(70)-mp.j) < mp.mpf('1e-65')
    values = [F(2,3),F(-4,5),F(7,11),F(3,2),F(1,7)]
    spectrum = finite_fourier(values,5)
    ctx = mp.mp.clone(); ctx.dps = 90
    for k,value in enumerate(spectrum):
        oracle = ctx.fsum(ctx.mpf(v.numerator)/v.denominator*ctx.exp(-2*ctx.pi*ctx.j*j*k/5)
                          for j,v in enumerate(values))/5
        assert abs(value.complex(80)-oracle) < ctx.mpf('1e-75')
    assert mp.mp.dps == before


def test_domains_resource_limits_and_unknown_cosets():
    K = Qp(3,10)
    f = TF(K,[1,2,3])
    assert f(K.O(1)) == 1
    with pytest.raises(PrecisionError): f(K.O(0))
    with pytest.raises(TypeError): f(Qp(5)(1))
    with pytest.raises(DomainError): f.translate(F(1,3))
    with pytest.raises(DomainError): f.convolve(TF.ball(K))
    with pytest.raises(TypeError): TF(K,[1.0,2,3])
    with pytest.raises(ValueError): TF(K,[1,2])
    with pytest.raises(ValueError): TF(K,[1],support_exponent=2,constancy_exponent=1)
    with pytest.raises(TypeError): TF(K,[1],support_exponent=True)
    with pytest.raises(PrecisionError): TF.ball(K,1,5,max_size=30)
    with pytest.raises(PrecisionError): finite_fourier(range(9),K,max_size=8)
    with pytest.raises(ValueError): finite_fourier([1,2],K)
    with pytest.raises(ValueError): finite_fourier([],K)
    with pytest.raises(TypeError): finite_fourier([1],K,inverse=1)
    with pytest.raises(ValueError): finite_fourier([1],K,max_size=True)
    with pytest.raises(PrecisionError): f.fourier(max_size=2)
    with pytest.raises(PrecisionError): f.convolve(f,max_size=2)
    with pytest.raises(TypeError): CV.root(3,1)+CV.root(5,1)
    with pytest.raises(ValueError): CV.root(3,1)**-1
    with pytest.raises(ZeroDivisionError): CV(3,[1])/0
    with pytest.raises(PrecisionError): CV.root(2,13)
    # Large-conductor phases remain available independently of root expansion.
    assert additive_phase(F(1,2**100),2) == F(1,2**100)
    with pytest.raises(PrecisionError): chi(F(1,2**100),2)
