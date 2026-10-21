"""p-adic analytic functions with checked domains and certified series tails.

The native algorithms prioritize transparent error bounds over asymptotic speed.
No series is stopped merely because one term happens to vanish modulo p**N.
"""
from __future__ import annotations

from fractions import Fraction

from .core import DomainError, Padic, PrecisionError, exact, p_power


def _target(x, precision):
    n = x.precision if precision is None else precision
    if type(n) is not int:
        raise TypeError("precision must be an integer")
    if n > x.precision:
        raise PrecisionError("function evaluation cannot invent input digits")
    return n


def _exp_domain(x):
    if not isinstance(x, Padic):
        raise TypeError("first construct a p-adic argument with Qp(p)(value)")
    if x.v * (x.p - 1) <= 1:
        raise DomainError("requires v_p(x) > 1/(p-1), including the whole input ball")


def _factorial_series(x, kind, precision=None):
    _exp_domain(x)
    nout = _target(x, precision)
    centre, term, total, k = x.representative(), Fraction(1), Fraction(0), 0
    # v(x^k/k!) >= k*(v(x)-1/(p-1)); a bound for the entire tail.
    slope = x.v * (x.p - 1) - 1
    while k == 0 or k * slope < nout * (x.p - 1):
        if kind == "exp":
            total += term
        elif kind == "sin" and k % 2:
            total += (-1)**((k-1)//2) * term
        elif kind == "cos" and k % 2 == 0:
            total += (-1)**(k//2) * term
        k += 1
        term = term * centre / k
    return x.field(total, absprec=nout)


def exp(x: Padic, *, precision=None) -> Padic:
    """Exponential on v_p(x)>1/(p-1); output absolute precision <= input."""
    return _factorial_series(x, "exp", precision)


def sin(x: Padic, *, precision=None) -> Padic:
    """Sine on the same sufficient convergence domain as exp."""
    return _factorial_series(x, "sin", precision)


def cos(x: Padic, *, precision=None) -> Padic:
    """Cosine on the same sufficient convergence domain as exp."""
    return _factorial_series(x, "cos", precision)


def _floor_log(n, p):
    r = 0
    while n >= p:
        n //= p
        r += 1
    return r


def log1p(x: Padic, *, precision=None) -> Padic:
    """log(1+x) on v_p(x)>=1, including p=2.

    The bound k*v_p(x)-floor(log_p(k)) is monotone for k>=1 and
    bounds every later term. This handles denominator spikes at p-powers.
    """
    if x.v < 1:
        raise DomainError("log1p requires v_p(x) >= 1")
    nout = _target(x, precision)
    c, power_x, total, k = x.representative(), x.representative(), Fraction(0), 1
    while k*x.v - _floor_log(k, x.p) < nout:
        total += (-1)**(k+1) * power_x / k
        k += 1
        power_x *= c
    return x.field(total, absprec=nout)


def teichmuller(x: Padic) -> Padic:
    """Root-of-unity part of a known unit.

    For odd p this lifts the residue to a (p-1)-st root of unity. At p=2
    returns the sign determined by the residue modulo 4; two digits are needed.
    """
    if x.is_zero() or x.v != 0:
        raise DomainError("teichmuller requires a known p-adic unit")
    n, p = x.precision, x.p
    if p == 2:
        if n < 2:
            raise PrecisionError("2-adic sign requires a residue modulo 4")
        return x.field(1 if x.residue(2) == 1 else -1, absprec=n)
    modulus, a = p**n, x.residue(1)
    for _ in range(n):
        a = pow(a, p, modulus)
    return x.field(a, absprec=n)


def log(x: Padic, *, p_branch=0, precision=None) -> Padic:
    """Iwasawa logarithm on Q_p^*, with the explicit branch log(p)=p_branch.

    Roots of unity have logarithm zero. The default branch is log(p)=0;
    log and exp are mutual inverses only on their common local domain.
    """
    if x.is_zero():
        raise PrecisionError("logarithm of a ball containing zero is undefined")
    v = x.v
    unit = x / p_power(x.p, v)
    principal = unit / teichmuller(unit)
    answer = log1p(principal - 1) + v*exact(p_branch)
    return answer if precision is None else answer.truncate(precision)


def power(base: Padic, exponent) -> Padic:
    """Principal-unit power exp(exponent*log(base)) for exponent in Z_p.

    Requires base in 1+p Z_p (odd p), or 1+4 Z_2. No implicit choice of
    fractional powers of roots of unity is made.
    """
    delta = base - 1
    _exp_domain(delta)
    if isinstance(exponent, Padic):
        base._check_prime(exponent)
        if exponent.v < 0:
            raise DomainError("exponent must lie in Z_p")
    else:
        from .core import valuation
        if valuation(exact(exponent), base.p) < 0:
            raise DomainError("exponent must lie in Z_p")
    return exp(log1p(delta) * exponent)


def _tonelli_shanks(a, p):
    if a % p == 0:
        return 0
    if pow(a, (p-1)//2, p) != 1:
        raise DomainError("unit is not a square modulo p")
    if p % 4 == 3:
        r = pow(a, (p+1)//4, p)
        return min(r, p-r)
    q, s = p-1, 0
    while q % 2 == 0:
        q //= 2
        s += 1
    z = 2
    while pow(z, (p-1)//2, p) != p-1:
        z += 1
    c, r, t, m = pow(z, q, p), pow(a, (q+1)//2, p), pow(a, q, p), s
    while t != 1:
        i, tt = 0, t
        while tt != 1:
            tt = tt*tt % p
            i += 1
        b = pow(c, 2**(m-i-1), p)
        r, t, c, m = r*b % p, t*b*b % p, b*b % p, i
    return min(r, p-r)


def sqrt(x: Padic) -> Padic:
    """One square root in Q_p, with the derivative's precision loss.

    Odd p: choose the smaller nonzero initial residue. At p=2 choose the
    unit root congruent to 1 modulo 4; a square unit must be 1 modulo 8.
    The negative gives the other root. Zero balls need further information.
    """
    if x.is_zero():
        raise PrecisionError("sqrt of a zero ball cannot select a root branch")
    if x.v % 2:
        raise DomainError("a square in Q_p has even valuation")
    k, rel, p, u = x.v//2, x.relative_precision, x.p, x.unit
    if p == 2:
        if rel < 3:
            raise PrecisionError("a 2-adic square unit needs three known digits")
        if u % 8 != 1:
            raise DomainError("a 2-adic square unit must be 1 modulo 8")
        r = 1
        for j in range(3, rel):
            r += (((u-r*r) >> j) & 1) * 2**(j-1)
        nout = x.precision - k - 1
    else:
        r, modulus = _tonelli_shanks(u % p, p), p
        for _ in range(1, rel):
            digit = ((u-r*r)//modulus * pow(2*r, -1, p)) % p
            r += digit * modulus
            modulus *= p
        nout = x.precision - k
    return x.field(r*p_power(p, k), absprec=nout)
