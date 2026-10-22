"""Exact coefficients for Fourier sums of prime-power roots of unity."""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from .core import DomainError, PrecisionError, exact, is_prime


def _size(p, level, limit):
    if not is_prime(p):
        raise ValueError("p must be prime")
    if type(level) is not int or level < 0:
        raise ValueError("level must be a nonnegative integer")
    if type(limit) is not int or limit < 1:
        raise ValueError("size limit must be a positive integer")
    size = 1
    for _ in range(level):
        size *= p
        if size > limit:
            raise PrecisionError(f"prime-power size exceeds the limit {limit}")
    return size


@dataclass(frozen=True, init=False, eq=False)
class CyclotomicValue:
    """An exact rational combination of p-power roots in the complex place.

    ``CyclotomicValue(p, coefficients, level=n)`` means sum(c_j*zeta**j),
    with zeta=exp(2*pi*i/p**n). Coefficients are rational, constant first.
    Reduction by Phi_(p**n), then descent to the smallest level, gives a
    canonical representation. Native root orders are limited to 4096.

    This is the exact coefficient ring for the harmonic APIs. It supports
    ring operations, conjugation and rational division; a general algebraic
    inverse or a p-adic extension-field embedding is not supplied here.
    """
    p: int
    level: int
    coefficients: tuple[Fraction, ...]

    def __init__(self, p, coefficients=(), *, level=0):
        result = self._from_terms(p,level,{
            i: exact(c) for i,c in enumerate(coefficients)})
        for key in ("p", "level", "coefficients"):
            object.__setattr__(self,key,getattr(result,key))

    @classmethod
    def _from_terms(cls, p, level, terms):
        order = _size(p,level,4096)
        reduced = {}
        for i,c in terms.items():
            index = i % order
            reduced[index] = reduced.get(index,Fraction(0))+c
        if level:
            step = order//p
            degree = (p-1)*step
            for i in tuple(reduced):
                if i >= degree:
                    c = reduced.pop(i)
                    for j in range(p-1):
                        index = i-degree+j*step
                        reduced[index] = reduced.get(index,Fraction(0))-c
        reduced = {i:c for i,c in reduced.items() if c}
        while level:
            if not reduced or set(reduced) <= {0}:
                level = 0
            elif level > 1 and all(i % p == 0 for i in reduced):
                reduced = {i//p:c for i,c in reduced.items()}
                level -= 1
            else:
                break
        size = max(reduced,default=-1)+1
        obj = object.__new__(cls)
        for key,value in (("p",p),("level",level),
                          ("coefficients",tuple(reduced.get(i,Fraction(0)) for i in range(size)))):
            object.__setattr__(obj,key,value)
        return obj

    @classmethod
    def root(cls, p, level=1, exponent=1):
        """The chosen exp(2*pi*i*exponent/p**level), represented exactly."""
        if type(exponent) is not int:
            raise TypeError("root exponent must be an integer")
        return cls._from_terms(p,level,{exponent:Fraction(1)})

    @property
    def order(self):
        """The root order of the canonical coefficient representation."""
        return self.p**self.level

    def as_fraction(self):
        """Return an exact rational value, raising if it is nonrational."""
        if self.level:
            raise DomainError("cyclotomic value is not rational")
        return self.coefficients[0] if self.coefficients else Fraction(0)

    def _coerce(self, other):
        if isinstance(other,CyclotomicValue):
            if other.p == self.p:
                return other
            if not other.level:
                return CyclotomicValue(self.p,[other.as_fraction()])
            raise TypeError("cyclotomic arithmetic requires the same prime")
        return CyclotomicValue(self.p,[exact(other)])

    def _terms(self, level):
        step = self.p**(level-self.level)
        return {i*step:c for i,c in enumerate(self.coefficients) if c}

    def __add__(self, other):
        other = self._coerce(other)
        level = max(self.level,other.level)
        terms = self._terms(level)
        for i,c in other._terms(level).items():
            terms[i] = terms.get(i,Fraction(0))+c
        return self._from_terms(self.p,level,terms)

    __radd__ = __add__

    def __neg__(self):
        return self._from_terms(self.p,self.level,{
            i:-c for i,c in enumerate(self.coefficients) if c})

    def __sub__(self, other):
        return self+-self._coerce(other)

    def __rsub__(self, other):
        return -self+other

    def __mul__(self, other):
        if not isinstance(other,CyclotomicValue):
            try:
                exact(other)
            except TypeError:
                return NotImplemented
        other = self._coerce(other)
        level = max(self.level,other.level)
        terms = {}
        for i,a in self._terms(level).items():
            for j,b in other._terms(level).items():
                terms[i+j] = terms.get(i+j,Fraction(0))+a*b
        return self._from_terms(self.p,level,terms)

    __rmul__ = __mul__

    def __truediv__(self, other):
        if isinstance(other,CyclotomicValue):
            other = other.as_fraction()
        return self*(1/exact(other))

    def __pow__(self, exponent):
        if type(exponent) is not int or exponent < 0:
            raise ValueError("cyclotomic exponent must be a nonnegative integer")
        result,base = CyclotomicValue(self.p,[1]),self
        while exponent:
            if exponent & 1:
                result *= base
            exponent >>= 1
            if exponent:
                base *= base
        return result

    def conjugate(self):
        """Exact complex conjugation for the chosen root of unity."""
        return self._from_terms(self.p,self.level,{
            -i:c for i,c in enumerate(self.coefficients) if c})

    def __eq__(self, other):
        if isinstance(other,CyclotomicValue):
            if not self.level and not other.level:
                return self.coefficients == other.coefficients
            return (self.p,self.level,self.coefficients) == (other.p,other.level,other.coefficients)
        try:
            return not self.level and self.as_fraction() == exact(other)
        except (TypeError,ValueError,ZeroDivisionError):
            return NotImplemented

    def __hash__(self):
        if not self.level:
            return hash(self.as_fraction())
        return hash((self.p,self.level,self.coefficients))

    def __bool__(self):
        return bool(self.coefficients)

    def complex(self, dps=50):
        """Approximate the chosen complex value (optional mpmath extra).

        dps controls working precision; this is not a certified error ball.
        The global mpmath context is unchanged.
        """
        from .bridges import _complex_context, _mp_rational
        ctx = _complex_context(dps)
        terms = [_mp_rational(ctx,c)*ctx.expjpi(ctx.mpf(2*i)/self.order)
                 for i,c in enumerate(self.coefficients) if c]
        return ctx.mpc(ctx.fsum(terms))

    def to_sympy(self):
        """An optional exact SymPy expression in the same complex place."""
        import sympy as sp
        root = sp.exp(2*sp.pi*sp.I/self.order)
        return sum((sp.Rational(c.numerator,c.denominator)*root**i
                    for i,c in enumerate(self.coefficients) if c),sp.S.Zero)
