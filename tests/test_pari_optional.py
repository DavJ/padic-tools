"""Independent optional PARI oracles for zeta, transforms and Taylor values.

Run where PARI/GP's gp command is installed. The normal exact-oracle tests
do not depend on it; skipping this file must never be reported as a pass
of the PARI cross-check.
"""
from fractions import Fraction
import shutil
import subprocess
import pytest
from padic_tools import AtomicMeasure, Qp, RationalAmiceMeasure, exp, finite_fourier, log, taylor, zeta_p


@pytest.mark.skipif(shutil.which("gp") is None, reason="optional PARI/GP executable not installed")
@pytest.mark.parametrize("p,s", [(3,-1),(5,-3),(7,-5),(3,2),(5,3),(7,4),(2,3)])
def test_principal_zeta_against_pari(p,s):
    script=f"print(lift(zeta({s}+O({p}^24))))\nquit\n"
    result=subprocess.run(["gp","-q","-f"],input=script,text=True,capture_output=True,check=True,timeout=20)
    oracle=Fraction(result.stdout.strip())
    assert zeta_p(s,Qp(p,12),branch=0).contains(oracle)


@pytest.mark.skipif(shutil.which("gp") is None, reason="optional PARI/GP executable not installed")
@pytest.mark.parametrize("p", [2,3,5])
def test_transforms_and_taylor_against_pari(p):
    K = Qp(p,12)
    point = 3 if p == 2 else 2
    denominator = 3 if p == 2 else 2
    omega = "-1" if p == 2 else "teichmuller(x)"
    script = (
        f"x={point}+O({p}^32); w={omega}; s=1/{denominator};\n"
        f"print(lift(exp(-{p*p}+O({p}^32))/{p}+2*exp(-2*{p*p}+O({p}^32))));\n"
        f"print(lift(w*exp(s*log(x/w))/{p}));\n"
        f"print(lift(exp({p*p+p**3}+O({p}^32))));\n"
        f"print(lift(log({p*p+p**3}+O({p}^32))));\nquit\n"
    )
    result = subprocess.run(["gp","-q","-f"],input=script,text=True,
                            capture_output=True,check=True,timeout=20)
    assert not result.stderr, result.stderr
    oracle = [Fraction(line) for line in result.stdout.splitlines()]
    assert len(oracle) == 4
    values = [
        AtomicMeasure(K,[(1,Fraction(1,p)),(2,2)]).laplace(p*p),
        AtomicMeasure(K,[(point,Fraction(1,p))]).mellin(Fraction(1,denominator),branch=1),
        taylor(exp,K(p*p),order=30)(K(p*p+p**3)),
        taylor(log,K(p*p),order=30)(K(p*p+p**3)),
    ]
    for value,reference in zip(values,oracle):
        assert value.precision >= 11
        assert value.contains(reference)


@pytest.mark.skipif(shutil.which("gp") is None, reason="optional PARI/GP executable not installed")
@pytest.mark.parametrize("p,c", [(2,3),(3,2),(5,2)])
def test_rational_measure_laplace_and_moments_against_pari(p,c):
    mu = RationalAmiceMeasure(Qp(p,14),[Fraction(1,p)],[c,1])
    script = (
        f"print(lift(1/({p}*({c}+exp(-{p*p}+O({p}^32))-1))));\n"
        f"f=1/({p}*({c}+exp(-x+O(x^10))-1));\n"
        "for(n=0,8,print(polcoef(f,n)))\nquit\n"
    )
    result = subprocess.run(["gp","-q","-f"],input=script,text=True,
                            capture_output=True,check=True,timeout=20)
    assert not result.stderr, result.stderr
    oracle = [Fraction(line) for line in result.stdout.splitlines()]
    assert len(oracle) == 10
    value = mu.laplace(p*p)
    assert value.precision >= 14
    assert value.contains(oracle[0])
    assert mu.laplace_series(9).coefficients == tuple(oracle[1:])
    from math import factorial
    for n,c in enumerate(oracle[1:]):
        assert mu.moment(n) == (-1)**n*factorial(n)*c


@pytest.mark.skipif(shutil.which("gp") is None, reason="optional PARI/GP executable not installed")
@pytest.mark.parametrize("p", [3,5,7])
def test_nonatomic_mellin_against_pari_kl_identity(p):
    # A(T)=1/(2+T), restricted to units, has the following analytic Mellin
    # identity on branch p-2. PARI evaluates its independent KL algorithm.
    mu = RationalAmiceMeasure(Qp(p,14),1,[2,1]).restrict_units()
    script = (
        f"s=1/2+O({p}^32);\n"
        f"print(lift((exp((s+1)*log(2+O({p}^32)))-1)*zeta(-s)));\nquit\n"
    )
    result = subprocess.run(["gp","-q","-f"],input=script,text=True,
                            capture_output=True,check=True,timeout=20)
    assert not result.stderr, result.stderr
    value = mu.mellin(Fraction(1,2),branch=p-2)
    assert value.precision == 14
    assert value.contains(Fraction(result.stdout.strip()))


@pytest.mark.skipif(shutil.which("gp") is None, reason="optional PARI/GP executable not installed")
@pytest.mark.parametrize("p,n", [(2,4),(3,9),(5,5)])
def test_exact_fourier_complex_place_against_pari(p,n):
    mp = pytest.importorskip("mpmath")
    script = (
        "default(realprecision,90);\n"
        f"N={n};\n"
        "for(k=0,N-1,z=sum(j=0,N-1,(j^2-2*j+3)*exp(-2*Pi*I*j*k/N))/N;"
        "print(real(z));print(imag(z)))\nquit\n"
    )
    result = subprocess.run(["gp","-q","-f"],input=script,text=True,
                            capture_output=True,check=True,timeout=20)
    assert not result.stderr, result.stderr
    oracle = result.stdout.splitlines()
    assert len(oracle) == 2*n
    ctx = mp.mp.clone(); ctx.dps = 85
    values = finite_fourier([j*j-2*j+3 for j in range(n)],p)
    for k,value in enumerate(values):
        reference = ctx.mpc(ctx.mpf(oracle[2*k].replace(" ","")),
                            ctx.mpf(oracle[2*k+1].replace(" ","")))
        assert abs(value.complex(80)-reference) < ctx.mpf("1e-74")
