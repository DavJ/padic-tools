"""Theta q-series and branches of the Kubota--Leopoldt p-adic zeta.

Conventions and the full tail argument are in docs/functions.rst.
"""
from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
from math import comb, gcd

from .core import DomainError, Padic, PadicField, PrecisionError, exact, valuation
from .analytic import exp, log, teichmuller


@lru_cache(maxsize=512)
def bernoulli(n: int) -> Fraction:
    """Exact Bernoulli number, with B_1=-1/2 (generating function t/(e^t-1))."""
    if type(n) is not int or n < 0:
        raise ValueError("n must be a nonnegative integer")
    if n == 0:
        return Fraction(1)
    if n > 1 and n % 2:
        return Fraction(0)
    return -sum((comb(n+1, k)*bernoulli(k) for k in range(n)), Fraction(0))/(n+1)


def zeta_special(s: int) -> Fraction:
    """Exact classical Riemann zeta at nonpositive integer arguments.

    This is algebraic special-value data, not a p-adic Dirichlet series.
    """
    if type(s) is not int or s > 0:
        raise DomainError("native exact special values require integer s <= 0")
    m = -s
    return (-1)**m * bernoulli(m+1)/(m+1)


def zeta_euler_value(s: int, p: int) -> Fraction:
    """Exact (1-p**(-s))*zeta(s) at integer s<=0, before selecting a branch."""
    PadicField(p)  # Validate p.
    return (1 - p**(-s))*zeta_special(s)


def zeta_p(s, field: PadicField | None = None, *, branch=0, precision=None) -> Padic:
    r"""Evaluate an even Kubota--Leopoldt branch on Z_p, away from s=1.

    Branch r interpolates (1-p**(n-1))*zeta(1-n) for n>=2 with
    n congruent to r modulo p-1 (p odd). At p=2 only r=0 is accepted,
    interpolating even n. r=0 is the usual principal zeta branch.

    The evaluator uses a Bernoulli expansion on residue discs, not a
    divergent p-adic sum of n**(-s). Finite input balls lose precision
    near the pole. A requested output precision is checked, never padded.
    Nonprincipal branches' removable value at s=1 is not implemented.
    """
    if isinstance(s, Padic):
        if field is not None and field.p != s.p:
            raise TypeError("field and argument have different primes")
        field = s.field if field is None else field
    elif field is None:
        raise TypeError("exact s requires a field, e.g. zeta_p(-3, Qp(5))")
    p, e = field.p, (2 if field.p == 2 else 1)
    period = 2 if p == 2 else p-1
    if type(branch) is not int or branch < 0 or branch >= period or branch % 2:
        raise DomainError("branch must be an even residue; at p=2 use branch=0")
    if isinstance(s, Padic):
        if s.v < 0:
            raise DomainError("s must lie in Z_p")
        delta = s-1
        if delta.is_zero():
            raise PrecisionError("input ball contains s=1; the pole cannot be resolved")
        d = delta.v
        target = min(field.precision, s.precision-e-2*d-2) if precision is None else precision
    else:
        q = exact(s)
        if valuation(q, p) < 0:
            raise DomainError("s must lie in Z_p")
        if q == 1:
            raise DomainError("s=1 is a pole on branch 0; other branch limits are unsupported")
        d = valuation(q-1, p)
        target = field.precision if precision is None else precision
    if type(target) is not int or target < 1:
        raise PrecisionError("at least one absolute digit is required; supply more input precision")
    work = target + e + 2*d + 5
    W = PadicField(p, work)
    arg = W(s) if isinstance(s, Padic) else W(q, absprec=work)
    F = p**e
    # v(B_k)>=-1 and binomial(1-s,k) in Z_p. After /F/(s-1),
    # every omitted term k>=K has valuation >= K*e - 1 - e - d.
    K = max(2, (target+1+e+d+e-1)//e)
    b, binomials = W(1, absprec=work), []
    for k in range(K):
        binomials.append(b)
        b = b * (1-arg-k)/(k+1)
    total = None
    for a in range(1, F):
        if gcd(a, p) != 1:
            continue
        aa = W(a, absprec=work)
        omega = teichmuller(aa)
        weight = exp((1-arg)*log(aa/omega)) * (omega**branch)
        inner = None
        for k, bk in enumerate(binomials):
            B = bernoulli(k)
            if not B:
                continue
            term = bk * (Fraction(F, a)**k * B)
            inner = term if inner is None else inner + term
        part = weight * inner
        total = part if total is None else total + part
    result = total / F / (arg-1)
    # Include the tail explicitly, even if cancellation in the partial sum
    # produces more digits than the proven bound.
    return result.truncate(target)


def theta3(q: Padic, z=1, *, precision=None) -> Padic:
    r"""Multiplicative theta sum: sum(n in Z) q**(n*n) * z**n.

    Requires v_p(q)>0 and nonzero z. For z=1 this is the classical theta
    constant 1+2*q+2*q**4+... evaluated in Q_p. For general z it matches
    Jacobi theta_3(w,q) when z=exp(2*i*w) over C; no such exponential
    identification is assumed p-adically.
    """
    if q.v <= 0:
        raise DomainError("theta requires |q|_p < 1 on the whole input ball")
    if isinstance(z, Padic):
        q._check_prime(z)
        if z.is_zero():
            raise PrecisionError("z must exclude zero because negative powers occur")
        w, zi = z.v, z.inverse()
    else:
        z = exact(z)
        if not z:
            raise DomainError("z cannot be zero")
        w, zi = valuation(z, q.p), 1/z
    # Changing q or z can amplify input error when z is nonintegral. Bound
    # the minimum of v(q)*n^2-|v(z)|*n over positive integer n, rather than
    # assuming a theta value has the precision of its nome.
    vertex = max(1, abs(w)//(2*q.v))
    min_term = min(q.v*j*j-abs(w)*j for j in (1, vertex, vertex+1))
    available = q.precision-q.v+min_term
    if isinstance(z, Padic):
        available = min(available, z.precision-w+min_term)
    target = min(q.field.precision, available) if precision is None else precision
    if type(target) is not int or target < 1:
        raise ValueError("theta output precision must be positive")
    n, total = 1, None
    while q.v*n*n - abs(w)*n < target or 2*q.v*n < abs(w):
        term = q**(n*n) * (z**n + zi**n)
        total = term if total is None else total + term
        n += 1
    if total is None:
        # The constant is exact; the uncertainty is the tail.
        return q.field(1, absprec=target)
    return (total+1).truncate(target)


def theta4(q: Padic, *, precision=None) -> Padic:
    """Theta constant 1+2*sum((-1)**n*q**(n*n), n>=1)."""
    return theta3(q, -1, precision=precision)
