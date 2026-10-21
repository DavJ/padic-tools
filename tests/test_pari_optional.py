"""Independent optional PARI oracle for the fixed principal branch.

Run where PARI/GP's gp command is installed. The normal exact-oracle tests
do not depend on it; skipping this file must never be reported as a pass
of the PARI cross-check.
"""
from fractions import Fraction
import shutil
import subprocess
import pytest
from padic_tools import Qp, zeta_p


@pytest.mark.skipif(shutil.which("gp") is None, reason="optional PARI/GP executable not installed")
@pytest.mark.parametrize("p,s", [(3,-1),(5,-3),(7,-5),(3,2),(5,3),(7,4),(2,3)])
def test_principal_zeta_against_pari(p,s):
    script=f"print(lift(zeta({s}+O({p}^24))))\nquit\n"
    result=subprocess.run(["gp","-q","-f"],input=script,text=True,capture_output=True,check=True,timeout=20)
    oracle=Fraction(result.stdout.strip())
    assert zeta_p(s,Qp(p,12),branch=0).contains(oracle)
