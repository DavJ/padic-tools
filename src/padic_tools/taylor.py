"""Exact formal Taylor jets and local p-adic Taylor enclosures."""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import comb, factorial

from . import analytic
from .core import DomainError, Padic, exact, valuation
from .symbolic import Polynomial, PowerSeries, exp_series, log1p_series
from .transforms import _order


@dataclass(frozen=True, init=False)
class TaylorExpansion:
    """A local p-adic expansion, constructed by ``taylor``.

    Coefficients are in ascending powers of (x-center). Calling the object
    evaluates the chosen truncation and adds a certified enclosure of the
    entire omitted tail. ``tail_precision(x)`` reports that remainder's
    absolute p-adic order, separately from coefficient/input uncertainty.
    """
    center: Padic
    coefficients: tuple
    _kind: str
    _scale: Padic | None
    _remainder: tuple

    def __init__(self, *args, **kwargs):
        raise TypeError("construct a TaylorExpansion through taylor()")

    @property
    def order(self):
        """First omitted Taylor degree."""
        return len(self.coefficients)

    def tail_precision(self, x):
        """Return B with remainder in O(p**B), or None for an exact polynomial.

        Raises outside the certified local convergence domain. A short
        expansion can legitimately give a coarse answer.
        """
        h = self.center.field(x) - self.center
        if self._kind == "polynomial":
            if not self._remainder:
                return None
            return min((c.v if isinstance(c, Padic) else valuation(c, h.p))
                       + n*h.v for n, c in enumerate(self._remainder, self.order))
        if self._kind in ("exp", "sin", "cos"):
            analytic._exp_domain(h)
            # For all k>=order: v(h**k/k!) >= k*(v(h)-1/(p-1)).
            slope = h.v*(h.p-1)-1
            return (self.order*slope + h.p-2) // (h.p-1)
        ratio = h / self._scale
        if ratio.v < 1:
            raise DomainError("log Taylor expansion requires v_p((x-center)/scale) >= 1")
        # k*v(ratio)-floor(log_p(k)) is nondecreasing for v(ratio)>=1.
        return self.order*ratio.v - analytic._floor_log(self.order, h.p)

    def __call__(self, x) -> Padic:
        """Evaluate with coefficient, input and analytic-tail uncertainty."""
        bound = self.tail_precision(x)
        h = self.center.field(x) - self.center
        value = self.coefficients[-1]
        for c in reversed(self.coefficients[:-1]):
            value = value*h + c
        if not isinstance(value, Padic):
            value = self.center.field(value)
        return value if bound is None else value + self.center.field.O(bound)


def taylor(function, center=0, *, order=10):
    """Taylor coefficients f^(n)(center)/n!, n=0,...,order-1.

    A Polynomial at an exact rational center returns a formal PowerSeries.
    At a Padic center, polynomials and exp/log/log1p/sin/cos return a
    TaylorExpansion with numerical tail bounds. The logarithm uses log(p)=0.
    Rational-centered standard exp/sin/cos/log1p jets at zero, and log at one,
    are also exact PowerSeries. Arbitrary numerical callables are rejected;
    derivatives are never inferred from finite differences.
    """
    _order(order)
    if not isinstance(center, Padic):
        a = exact(center)
        if isinstance(function, Polynomial):
            return PowerSeries(function.shift(a).coefficients, order)
        if function is analytic.log and a == 1:
            return log1p_series(order)
        if a == 0:
            if function is analytic.exp:
                return exp_series(order)
            if function is analytic.log1p:
                return log1p_series(order)
            if function is analytic.sin:
                return PowerSeries([Fraction((-1)**((n-1)//2), factorial(n))
                                    if n % 2 else 0 for n in range(order)], order)
            if function is analytic.cos:
                return PowerSeries([Fraction((-1)**(n//2), factorial(n))
                                    if n % 2 == 0 else 0 for n in range(order)], order)
        raise TypeError("this analytic center requires an explicit p-adic ball")

    scale, remainder = None, ()
    if isinstance(function, Polynomial):
        kind = "polynomial"
        cs = []
        for n in range(function.degree+1):
            cs.append(sum((function[j]*comb(j,n) *
                           (center**(j-n) if j != n else 1)
                           for j in range(n, function.degree+1)), Fraction(0)))
        remainder = tuple(cs[order:])
        cs = cs[:order] + [Fraction(0)]*max(0, order-len(cs))
    elif function in (analytic.exp, analytic.sin, analytic.cos):
        kind = function.__name__
        value = function(center)  # Also checks the whole center ball.
        if function is analytic.exp:
            derivatives = (value,)
        else:
            sine, cosine = analytic.sin(center), analytic.cos(center)
            derivatives = ((sine, cosine, -sine, -cosine) if function is analytic.sin
                           else (cosine, -sine, -cosine, sine))
        cs = [derivatives[n % len(derivatives)] / factorial(n) for n in range(order)]
    elif function in (analytic.log, analytic.log1p):
        kind = "log"
        scale = center if function is analytic.log else 1+center
        cs = [function(center)]
        cs.extend(((-1)**(n+1) / (scale**n)) / n for n in range(1,order))
    else:
        raise TypeError("supported functions are Polynomial, exp, log, log1p, sin and cos")
    obj = object.__new__(TaylorExpansion)
    for name, value in (("center", center), ("coefficients", tuple(cs)),
                        ("_kind", kind), ("_scale", scale), ("_remainder", remainder)):
        object.__setattr__(obj, name, value)
    return obj
