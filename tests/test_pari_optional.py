"""Independent optional PARI oracles for zeta, transforms and Taylor values.

Run where PARI/GP's gp command is installed. The normal exact-oracle tests
do not depend on it; skipping this file must never be reported as a pass
of the PARI cross-check.
"""
from fractions import Fraction
import shutil
import subprocess
import pytest
from padic_tools import AtomicMeasure, Qp, exp, log, taylor, zeta_p


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
