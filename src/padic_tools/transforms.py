"""p-adic transforms of finite rational atomic measures on Z_p.

The measure is explicit data, not an implicit Haar or Volkenborn integral.
Numerical kernels use ball arithmetic; formal series retain exact moments.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import factorial

from .analytic import exp, power, teichmuller
from .core import DomainError, Padic, PadicField, exact, valuation
from .symbolic import PowerSeries


def _order(order):
    if type(order) is not int or order < 1:
        raise ValueError("series order must be a positive integer")
    return order


@dataclass(frozen=True, init=False)
class AtomicMeasure:
    """A bounded measure ``sum(weight * delta(point))`` on Z_p.

    ``AtomicMeasure(field, [(point, weight), ...])`` accepts exact rational
    data. Points must be p-integral; weights may have negative valuation.
    Repeated points are combined and exact zero weights removed. Negative
    weights are allowed: this is not necessarily a probability measure.
    """
    field: PadicField
    atoms: tuple[tuple[Fraction, Fraction], ...]

    def __init__(self, field: PadicField, atoms=()):
        if not isinstance(field, PadicField):
            raise TypeError("field must be Qp(p, precision)")
        combined = {}
        for point, weight in atoms:
            point, weight = exact(point), exact(weight)
            if valuation(point, field.p) < 0:
                raise DomainError("atomic measure points must lie in Z_p")
            combined[point] = combined.get(point, Fraction(0)) + weight
        object.__setattr__(self, "field", field)
        object.__setattr__(self, "atoms", tuple(sorted(
            (point, weight) for point, weight in combined.items() if weight)))

    def moment(self, n: int) -> Fraction:
        """The exact moment integral x**n dmu, for a nonnegative integer n."""
        if type(n) is not int or n < 0:
            raise ValueError("moment index must be a nonnegative integer")
        return sum((weight * point**n for point, weight in self.atoms), Fraction(0))

    def restrict_units(self) -> AtomicMeasure:
        """Explicit restriction to Z_p^*, discarding all nonunit atoms."""
        return AtomicMeasure(self.field, (
            (point, weight) for point, weight in self.atoms
            if valuation(point, self.field.p) == 0))

    def laplace_series(self, order=10) -> PowerSeries:
        """Exact formal series integral exp(-t*x) dmu modulo t**order.

        Coefficient n is (-1)**n * moment(n) / n!. This method alone makes
        no numerical remainder assertion; ``laplace`` checks convergence.
        """
        return PowerSeries([(-1)**n * self.moment(n) / factorial(n)
                            for n in range(_order(order))], order)

    def amice_series(self, order=10) -> PowerSeries:
        """Exact formal Amice transform integral (1+T)**x dmu.

        Coefficient n is the binomial moment integral binomial(x,n) dmu.
        """
        cs = [Fraction(0)] * _order(order)
        for point, weight in self.atoms:
            choose = Fraction(1)
            for n in range(order):
                cs[n] += weight * choose
                choose *= (point - n) / (n + 1)
        return PowerSeries(cs, order)

    def amice(self, parameter) -> Padic:
        """Evaluate the Amice transform on v_p(T)>=1, including p=2.

        Binomial coefficients of p-integral points are p-integral. If b is
        a lower bound for all weight valuations, the omitted tail from n
        has valuation at least b+n*v_p(T). Input precision is propagated
        separately through the polynomial evaluation.
        """
        t = self.field(parameter)
        if t.v < 1:
            raise DomainError("Amice evaluation requires v_p(T) >= 1")
        if not self.atoms:
            return self.field.O(t.precision)
        b = min(valuation(weight, t.p) for _, weight in self.atoms)
        order = max(1, (t.precision + t.v - 1) // t.v)
        value = self.amice_series(order).polynomial()(t)
        return value + self.field.O(b + order*t.v)

    def laplace(self, parameter) -> Padic:
        """Local Laplace transform integral exp(-t*x) dmu.

        Every nonzero atom must satisfy v_p(t*x)>1/(p-1). The whole input
        ball is checked. The minus sign is the chosen Laplace convention.
        """
        t = self.field(parameter)
        value = Fraction(0)
        for point, weight in self.atoms:
            value += weight * exp(-t * point) if point else weight
        return value if isinstance(value, Padic) else self.field(value, absprec=t.precision)

    def mellin(self, exponent, *, branch=0) -> Padic:
        """Mellin transform integral omega(x)**branch * <x>**s dmu.

        Support must lie in Z_p^*. Exponent s lies in Z_p; <x>=x/omega(x).
        Branch is modulo p-1 for odd p, modulo 2 for the 2-adic sign.
        At nonnegative integers n congruent to branch, this equals the exact n-th
        moment. Restriction of a measure to units is always explicit.
        """
        if type(branch) is not int:
            raise TypeError("Mellin branch must be an integer")
        s = self.field(exponent)
        if s.v < 0:
            raise DomainError("Mellin exponent must lie in Z_p")
        if any(valuation(point, s.p) != 0 for point, _ in self.atoms):
            raise DomainError("Mellin support must lie in Z_p^*; use restrict_units()")
        branch %= 2 if s.p == 2 else s.p - 1
        # Atoms are exact: guard digits can be requested without extending s.
        work_precision = max(self.field.precision, s.precision) + 2
        value = Fraction(0)
        for point, weight in self.atoms:
            x = self.field(point, absprec=work_precision)
            omega = teichmuller(x)
            torsion = omega**branch if branch else 1
            value += weight * torsion * power(x / omega, s)
        return value if isinstance(value, Padic) else self.field(value, absprec=s.precision)
