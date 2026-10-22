"""Bounded Z_p measures specified by exact rational Amice transforms.

Integral normalized denominator coefficients prove a uniform bound on the
binomial moments. No finite atom list or sampling cutoff is needed.
"""
from __future__ import annotations

from dataclasses import dataclass, field as data_field
from fractions import Fraction
from functools import lru_cache
from math import factorial

from .analytic import exp
from .core import DomainError, Padic, PadicField, PrecisionError, exact, valuation
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
    assuming that a finite prefix has zero tail. Unit restriction is exact;
    Mellin evaluation on unit-supported measures uses certified Mahler tails.
    """
    field: PadicField
    numerator: Polynomial
    denominator: Polynomial
    valuation_bound: int
    _unit_support: bool = data_field(default=False,repr=False,compare=False)

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
                           ("valuation_bound",b),("_unit_support",False)):
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
        result = cls(measure.field,numerator,unit**offset)
        if all(valuation(point,measure.field.p) == 0 for point,_ in measure.atoms):
            object.__setattr__(result,"_unit_support",True)
        return result

    def restrict_units(self, *, max_degree=4096):
        """Exact restriction to Z_p^*, returning a rational Amice measure.

        In z=1+T, subtract the Laurent-series section with exponents
        divisible by p. Decimated recurrences compute this rational function
        over Q, without numerical p-th roots of unity. max_degree limits
        intermediate polynomial size; the native recurrence order is <=64.
        """
        _order(max_degree)
        if self._unit_support or self.numerator.degree < 0:
            return self
        return self._restrict_units(max_degree)

    @lru_cache(maxsize=64)
    def _restrict_units(self, max_degree):
        from ._rational_sections import zero_section
        numerator,denominator = self.numerator.shift(-1),self.denominator.shift(-1)
        top,bottom = zero_section(numerator,denominator,self.field.p,max_degree)
        result = RationalAmiceMeasure(self.field,
            (numerator*bottom-top*denominator).shift(1),
            (denominator*bottom).shift(1))
        object.__setattr__(result,"_unit_support",True)
        return result

    def mellin(self, exponent, *, branch=0, precision=None, max_terms=10000):
        """Integral omega(x)**branch * <x>**s dmu on unit support.

        Restriction must be explicit if the measure has nonunit support.
        Absolute output precision defaults to the field precision, capped by
        existing exponent digits (with the derivative's gain). Exact exponent
        data may be evaluated at additional working precision. max_terms is a
        resource limit; insufficient digits or terms raise PrecisionError.

        Mahler coefficients of the kernel, extended by zero on nonunits,
        have valuation >=v_p(n!) at odd p and >=floor(n/2)-1 at p=2.
        Multiplication by the measure adds valuation_bound to the tail order.
        """
        if type(branch) is not int:
            raise TypeError("Mellin branch must be an integer")
        _order(max_terms)
        target = self.field.precision if precision is None else precision
        if type(target) is not int:
            raise TypeError("precision must be an integer")
        s = self.field(exponent)
        if s.v < 0:
            raise DomainError("Mellin exponent must lie in Z_p")
        if not self._unit_support and self.restrict_units() != self:
            raise DomainError("Mellin support must lie in Z_p^*; use restrict_units()")
        p,b = self.field.p,self.valuation_bound
        gain = 2 if p == 2 else 1
        if isinstance(exponent,Padic):
            available = s.precision+gain+b
            if precision is not None and target > available:
                raise PrecisionError("Mellin precision exceeds the available exponent digits")
            target = min(target,available)
        needed = target-b
        if needed <= 0 or self.numerator.degree < 0:
            return self.field.O(target)
        digits = max(0,needed-gain)
        if not isinstance(exponent,Padic) and s.precision < digits:
            s = self.field(exponent,absprec=digits)
        if p == 2:
            count = 2*(needed+1)
        else:
            lo,hi = 1,p*needed
            while lo < hi:
                mid = (lo+hi)//2
                n,v = mid,0
                while n:
                    n //= p
                    v += n
                if v >= needed:
                    hi = mid
                else:
                    lo = mid+1
            count = lo
        if count > max_terms:
            raise PrecisionError(f"Mellin needs {count} terms; max_terms={max_terms}")
        modulus = p**needed
        residue = s.residue(digits)
        branch %= 2 if p == 2 else p-1
        values = []
        if p != 2:
            # CRT chooses n=s mod p**digits and n=branch mod p-1.
            power = residue+p**digits*((branch-residue) % (p-1))
        for k in range(count):
            if k % p == 0:
                values.append(0)
            elif p == 2:
                sign = 1 if k % 4 == 1 else -1
                values.append(sign**branch*pow(sign*k,residue,modulus) % modulus)
            else:
                values.append(pow(k,power,modulus))
        # Forward differences compute Mahler coefficients; integer operations
        # preserve the uniform modulus even when exponent errors are correlated.
        moments = self.amice_series(count).coefficients
        value = Fraction(0)
        for moment in moments:
            if moment:
                value += self.field(values[0],absprec=needed)*moment
            values = [(right-left) % modulus for left,right in zip(values,values[1:])]
        if not isinstance(value,Padic):
            value = self.field(value,absprec=target)
        return value + self.field.O(target)

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
        result = RationalAmiceMeasure(self.field,
            self.numerator*other.denominator+other.numerator*self.denominator,
            self.denominator*other.denominator)
        if self._unit_support and other._unit_support:
            object.__setattr__(result,"_unit_support",True)
        return result

    def scale(self, weight):
        """Multiply all measure values by an exact rational weight."""
        result = RationalAmiceMeasure(self.field,self.numerator*exact(weight),self.denominator)
        if self._unit_support:
            object.__setattr__(result,"_unit_support",True)
        return result

    def convolve(self, other):
        """Additive convolution: Amice transforms multiply."""
        self._check(other)
        return RationalAmiceMeasure(self.field,self.numerator*other.numerator,
                                   self.denominator*other.denominator)
