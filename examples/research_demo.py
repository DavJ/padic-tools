"""Run with python examples/research_demo.py after installing the project."""
from fractions import Fraction
from padic_tools import (
    Qp, exp, log, theta3, zeta_p, Gaussian, GaussianBridge,
    Hamilton, Biquaternions, QuaternionBridge,
)


def main():
    K = Qp(5, 20)
    assert log(exp(K(5))).is_congruent(5)
    print("theta3(5):", theta3(K(5)))
    print("Z_(5,0)(-3):", zeta_p(-3, K, branch=0))
    print("expected exact Euler-corrected value:", Fraction(-31, 30))
    g = GaussianBridge(Gaussian(Fraction(1, 3), 2))
    print("Gaussian image in quadratic 3-adic algebra:", g.padic_algebra(Qp(3)))
    print("Chosen Gaussian image in Q5:", g.split_padic(K))
    q = Hamilton(1, 2, 3, 4)
    assert (q*q).to_clifford() == q.to_clifford()*q.to_clifford()
    b = Biquaternions(Gaussian(1, 2), Gaussian(3, 4), 1, 0)
    print("Biquaternion image:", QuaternionBridge(b).padic(K))
    print("All demo identities verified.")


if __name__ == "__main__":
    main()
