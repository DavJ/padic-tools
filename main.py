"""Compatibility helpers for the 2022 prototype; use padic_tools for new code."""
from padic_tools import Qp, exact


def rationalize(number):
    """Return exact numerator and denominator; binary floats are rejected."""
    value = exact(number)
    return value.numerator, value.denominator


def bezout2(a: int, b: int):
    """Return integers (s,t,g) with s*a+t*b=g=gcd(a,b), including zero inputs."""
    if type(a) is not int or type(b) is not int:
        raise TypeError("Bezout inputs must be integers")
    r0, r1, s0, s1, t0, t1 = a, b, 1, 0, 0, 1
    while r1:
        q = r0//r1
        r0, r1 = r1, r0-q*r1
        s0, s1 = s1, s0-q*s1
        t0, t1 = t1, t0-q*t1
    sign = -1 if r0 < 0 else 1
    return sign*s0, sign*t0, sign*r0


def bezout1(a: int, b: int):
    """Return the two integer Bezout coefficients."""
    return bezout2(a, b)[:2]


def padic(nominator, denominator=1, prime=5, precision=20):
    """The original unfinished entry point, now returning a Padic ball."""
    return Qp(prime, precision)(exact(nominator)/exact(denominator))


if __name__ == "__main__":
    print(padic(1, 3))
