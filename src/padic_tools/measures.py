"""Bounded Z_p measures specified by exact rational Amice transforms.

Integral normalized denominator coefficients prove a uniform bound on the
binomial moments. No finite atom list or sampling cutoff is needed.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from math import factorial

from .analytic import exp
from .core import DomainError, Padic, PadicField, exact, valuation
from .symbolic import Polynomial, PowerSeries
from .transforms import AtomicMeasure, _order


def _polynomial(value):
    if isinstance(value, Polynomial):
        return value
    if isinstance(value, (list, tuple)):
        return Polynomial(value)
    return Polynomial([value])


def _evaluate(polynomial, x):
    # Keep constant polynomials exact instead of rounding at field defaults.
    if not polynomial.coefficients:
        return Fraction(0)
    value = polynomial.coefficients[-1]
    for c in reversed(polynomial.coefficients[:-1]):
        value = value*x+c
    return value


def _index(n):
    if type(n) is not int or n < 0:
        raise ValueError("moment index must be a nonnegative integer")
    return n


@dataclass(frozen=True, init=False)
class RationalAmiceMeasure:
    """A bounded Qp-valued measure on Z_p with transform P(T)/Q(T).

    Numerator and denominator are exact rational polynomials, coefficient
    lists, or scalars. Common factors are cancelled, then Q is normalized
    to Q(0)=1. Every other coefficient of Q must be p-integral. This proves
    that every Amice coefficient has valuation at least ``valuation_bound``.

    This represents an infinite sequence of binomial moments, rather than
    assuming that a finite prefix has zero tail. General Mellin evaluation
    on units and inverse transforms are not part of this class yet.
    """
    field: PadicField
    numerator: Polynomial
    denominator: Polynomial
    valuation_bound: int

    def __init__(self, field, numerator, denominator=1):
        if not isinstance(field, PadicField):
            raise TypeError("field must be Qp(p, precision)")
        p,q = _polynomial(numerator),_polynomial(denominator)
        if q.degree < 0:
            raise ZeroDivisionError("zero Amice denominator")
        common = p.gcd(q)
        p,q = divmod(p,common)[0],divmod(q,common)[0]
        if not q[0]:
            raise DomainError("reduced Amice transform has a pole at T=0")
        p,q = p*(1/q[0]),q*(1/q[0])
        if any(valuation(c,field.p) < 0 for c in q.coefficients):
            raise DomainError("normalized Amice denominator coefficients must be p-integral")
        b = min((valuation(c,field.p) for c in p.coefficients if c), default=0)
        for name,value in (("field",field),("numerator",p),("denominator",q),
                           ("valuation_bound",b)):
            object.__setattr__(self,name,value)

    @classmethod
    def from_atomic(cls, measure):
        """Convert integer-supported AtomicMeasure data exactly.

        Negative integers give denominator powers of (1+T). Noninteger
        rational support generally gives algebraic, not rational, transforms
        and is rejected rather than approximated.
        """
        if not isinstance(measure, AtomicMeasure):
            raise TypeError("expected an AtomicMeasure")
        if any(point.denominator != 1 for point,_ in measure.atoms):
            raise DomainError("rational Amice conversion requires integer support")
        offset = max((max(0,-int(point)) for point,_ in measure.atoms), default=0)
        unit = Polynomial([1,1])
        numerator = Polynomial()
        for point,weight in measure.atoms:
            numerator += weight * unit**(int(point)+offset)
        return cls(measure.field,numerator,unit**offset)

    def amice_series(self, order=10):
        """Exact first ``order`` binomial moments as a formal PowerSeries."""
        return self._series(_order(order))

    @lru_cache(maxsize=128)
    def _series(self, order):
        p,q = self.numerator,self.denominator
        cs = []
        for n in range(order):
            cs.append(p[n]-sum((q[j]*cs[n-j] for j in range(1,min(n,q.degree)+1)),
                               Fraction(0)))
        return PowerSeries(cs,order)

    def binomial_moment(self, n):
        """Exact integral binomial(x,n) dmu."""
        return self.amice_series(_index(n)+1).coefficients[n]

    def moment(self, n):
        """Exact integral x**n dmu, using Stirling numbers of the second kind."""
        _index(n)
        row = [1]
        for degree in range(1,n+1):
            row = [0]+[row[k-1]+k*(row[k] if k < len(row) else 0)
                       for k in range(1,degree+1)]
        a = self.amice_series(n+1).coefficients
        return sum((row[k]*factorial(k)*a[k] for k in range(n+1)),Fraction(0))

    def integrate_polynomial(self, function):
        """Exact integral of a rational Polynomial on Z_p."""
        if not isinstance(function, Polynomial):
            raise TypeError("expected a Polynomial over Q")
        return sum((c*self.moment(n) for n,c in enumerate(function.coefficients)),Fraction(0))

    def integrate_mahler(self, coefficients, *, tail_valuation):
        """Pair an exact Mahler prefix with a user-supplied uniform tail bound.

        If f=sum c_n*binomial(x,n) is continuous and every omitted c_n has
        valuation >= B=tail_valuation, the omitted integral has valuation
        >= B+valuation_bound. Establishing B and continuity is the caller's
        responsibility: a finite prefix does not prove an infinite bound.
        """
        if type(tail_valuation) is not int:
            raise TypeError("tail_valuation must be an integer")
        cs = tuple(exact(c) for c in coefficients)
        moments = self.amice_series(max(1,len(cs))).coefficients
        value = sum((c*a for c,a in zip(cs,moments)),Fraction(0))
        return self.field(value,absprec=tail_valuation+self.valuation_bound)

    def amice(self, parameter):
        """Evaluate P(T)/Q(T) on v_p(T)>=1, propagating input precision.

        The normalization proves that Q(T) is a unit throughout this disc.
        Rational evaluation includes the whole series; there is no cutoff.
        """
        t = self.field(parameter)
        if t.v < 1:
            raise DomainError("Amice evaluation requires v_p(T) >= 1")
        value = _evaluate(self.numerator,t)/_evaluate(self.denominator,t)
        return value if isinstance(value,Padic) else self.field(
            value,absprec=t.precision+self.valuation_bound)

    def laplace(self, parameter):
        """Local integral exp(-t*x) dmu for v_p(t)>1/(p-1)."""
        t = self.field(parameter)
        return self.amice(exp(-t)-1)

    def laplace_series(self, order=10):
        """Exact formal Laplace coefficients, via A_mu(exp(-t)-1)."""
        t = PowerSeries([0,1],_order(order))
        return self.amice_series(order).compose((-t).exp()-1)

    def _check(self, other):
        if not isinstance(other,RationalAmiceMeasure):
            raise TypeError("expected a RationalAmiceMeasure")
        if other.field.p != self.field.p:
            raise TypeError("measure operations require the same prime")

    def __add__(self, other):
        if not isinstance(other,RationalAmiceMeasure):
            return NotImplemented
        self._check(other)
        return RationalAmiceMeasure(self.field,
            self.numerator*other.denominator+other.numerator*self.denominator,
            self.denominator*other.denominator)

    def scale(self, weight):
        """Multiply all measure values by an exact rational weight."""
        return RationalAmiceMeasure(self.field,self.numerator*exact(weight),self.denominator)

    def convolve(self, other):
        """Additive convolution: Amice transforms multiply."""
        self._check(other)
        return RationalAmiceMeasure(self.field,self.numerator*other.numerator,
                                   self.denominator*other.denominator)
