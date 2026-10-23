"""Independent optional PARI oracles for zeta, transforms and Taylor values.

Run where PARI/GP's gp command is installed. The normal exact-oracle tests
do not depend on it; skipping this file must never be reported as a pass
of the PARI cross-check.
"""
from fractions import Fraction
import shutil
import subprocess
import pytest
from padic_tools import AtomicMeasure, BallVladimirov, EisensteinQuadraticField, Gaussian, GaussianLocalPlace, PadicTestFunction, Qp, RationalAmiceMeasure, UnramifiedQuadraticField, exp, finite_fourier, log, taylor, zeta_p


@pytest.mark.skipif(shutil.which("gp") is None, reason="optional PARI/GP executable not installed")
@pytest.mark.parametrize("p,sign", [(5, 1), (5, -1), (13, 1), (13, -1)])
def test_split_gaussian_place_and_guard_digits_against_pari(p, sign):
    script = (
        f"p={p};r=sqrt(-1+O(p^60));if(lift(r)%p>p/2,r=-r);\n"
        f"z=1/p^6-2*{sign}*r/p^5;\n"
        f"print(lift(z));print(lift(1/z));print(lift({sign}*r));quit\n"
    )
    result = subprocess.run(["gp", "-q", "-f"], input=script, text=True,
                            capture_output=True, check=True, timeout=20)
    assert not result.stderr, result.stderr
    oracle = [Fraction(line) for line in result.stdout.splitlines()]
    assert len(oracle) == 3
    place = GaussianLocalPlace(Qp(p, 4), root_sign=sign)
    image = place(Gaussian(Fraction(1, p**6), Fraction(-2, p**5)), absprec=11)
    assert image.precision == 11 and image.contains(oracle[0])
    assert image.inverse().contains(oracle[1])
    assert place.generator_image(absprec=13).contains(oracle[2])


@pytest.mark.skipif(shutil.which("gp") is None, reason="optional PARI/GP executable not installed")
@pytest.mark.parametrize("p,sign", [(2, 1), (2, -1), (3, 1), (3, -1)])
def test_nonsplit_gaussian_places_against_pari_quotients(p, sign):
    polynomial, root = ("t^2-2*t+2", "(t-1)") if p == 2 else ("t^2+1", "t")
    script = (
        f"p={p};z=Mod((1/p^3+O(p^40))+{sign}*(2/p^2+O(p^40))*{root},{polynomial});\n"
        "print(lift(polcoef(lift(z),0)));print(lift(polcoef(lift(z),1)));\n"
        "w=1/z;print(lift(polcoef(lift(w),0)));print(lift(polcoef(lift(w),1)));\n"
        "print(lift(norm(z)));print(lift(trace(z)));quit\n"
    )
    result = subprocess.run(["gp", "-q", "-f"], input=script, text=True,
                            capture_output=True, check=True, timeout=20)
    assert not result.stderr, result.stderr
    oracle = [Fraction(line) for line in result.stdout.splitlines()]
    assert len(oracle) == 6
    place = GaussianLocalPlace(Qp(p, 4), root_sign=sign)
    image = place(Gaussian(Fraction(1, p**3), Fraction(2, p**2)), absprec=19)
    assert image.precision == 19 and image.contains(tuple(oracle[:2]))
    assert image.inverse().contains(tuple(oracle[2:4]))
    assert image.norm().contains(oracle[4]) and image.trace().contains(oracle[5])


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


@pytest.mark.skipif(shutil.which("gp") is None, reason="optional PARI/GP executable not installed")
@pytest.mark.parametrize("p,a,b,order", [(2,0,2,1),(3,-1,1,2),(5,0,1,1)])
def test_regional_vladimirov_against_independent_pari_kernel(p,a,b,order):
    n = p**(b-a)
    script = (
        f"p={p};a={a};b={b};s={order};N={n};C=(p^s-1)/(1-p^(-1-s));\n"
        "M=matrix(N,N,j,k,if(j==k,0,-C*p^(-b)*p^((1+s)*(a+valuation(j-k,p)))));\n"
        "for(j=1,N,M[j,j]=-sum(k=1,N,M[j,k]));\n"
        "v=vector(N,j,(j-1)^2+1)~;w=M*v;for(j=1,N,print(w[j]));\nquit\n"
    )
    result = subprocess.run(["gp","-q","-f"],input=script,text=True,
                            capture_output=True,check=True,timeout=20)
    assert not result.stderr, result.stderr
    oracle = tuple(Fraction(line) for line in result.stdout.splitlines())
    assert len(oracle) == n
    op = BallVladimirov(p,order,support_exponent=a,constancy_exponent=b)
    f = PadicTestFunction(p,[j*j+1 for j in range(n)],support_exponent=a,constancy_exponent=b)
    assert tuple(value.as_fraction() for value in op.apply(f).values) == oracle


@pytest.mark.skipif(shutil.which("gp") is None, reason="optional PARI/GP executable not installed")
@pytest.mark.parametrize("p,c,b", [(2,1,1),(3,1,0),(5,2,0),(7,1,0)])
def test_unramified_quadratic_against_pari_polymod(p,c,b):
    # PARI inverts a polynomial modulo f over its own Qp coefficient field;
    # it independently computes the norm and trace of that polymod.
    script = (
        f"f=t^2+{b}*t+{c};z=Mod((1/{p}+O({p}^32))+(3+O({p}^32))*t,f);\n"
        "w=lift(z^-1);print(lift(polcoef(w,0)));print(lift(polcoef(w,1)));\n"
        "print(lift(norm(z)));print(lift(trace(z)));\n"
        f"w=lift(z*Mod(1+{p}*t,f));print(lift(polcoef(w,0)));print(lift(polcoef(w,1)));\n"
        f"print(valuation(lift(z),{p}));quit\n"
    )
    result = subprocess.run(["gp","-q","-f"],input=script,text=True,
                            capture_output=True,check=True,timeout=20)
    assert not result.stderr, result.stderr
    oracle = [Fraction(line) for line in result.stdout.splitlines()]
    assert len(oracle) == 7
    L = UnramifiedQuadraticField(Qp(p,12),[c,b,1])
    x = L(Fraction(1,p),3)
    assert x.inverse().contains(oracle[:2])
    assert x.norm().contains(oracle[2])
    assert x.trace().contains(oracle[3])
    assert (x*L(1,p)).contains(oracle[4:6])
    assert x.valuation() == oracle[6]


@pytest.mark.skipif(shutil.which("gp") is None, reason="optional PARI/GP executable not installed")
@pytest.mark.parametrize("p,c,b", [(2,2,-2),(2,-2,0),(3,-3,0),(5,5,5)])
def test_eisenstein_quadratic_against_pari_polymod_and_norm_valuation(p,c,b):
    points = [(0,1),(Fraction(1,p),3),(0,Fraction(1,p)),(p,p)]
    script = f"f=t^2+({b})*t+({c});\n"
    for a,d in points:
        script += (
            f"z=Mod(({a}+O({p}^40))+({d}+O({p}^40))*t,f);\n"
            "w=lift(z^-1);print(lift(polcoef(w,0)));print(lift(polcoef(w,1)));\n"
            "print(lift(norm(z)));print(lift(trace(z)));\n"
            f"w=lift(z*Mod(1+{p}*t,f));print(lift(polcoef(w,0)));print(lift(polcoef(w,1)));\n"
            f"print(valuation(norm(z),{p}));\n"
        )
    result = subprocess.run(["gp","-q","-f"],input=script+"quit\n",text=True,
                            capture_output=True,check=True,timeout=20)
    assert not result.stderr, result.stderr
    oracle = [Fraction(line) for line in result.stdout.splitlines()]
    assert len(oracle) == 7*len(points)
    L = EisensteinQuadraticField(Qp(p,12),[c,b,1])
    for j,(a,d) in enumerate(points):
        x, reference = L(a,d), oracle[7*j:7*j+7]
        assert x.inverse().contains(reference[:2])
        assert x.norm().contains(reference[2])
        assert x.trace().contains(reference[3])
        assert (x*L(1,p)).contains(reference[4:6])
        assert x.valuation() == reference[6]
