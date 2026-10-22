"""Exact finite harmonic analysis of locally constant functions on Q_p.

The coefficient domain is rational combinations of complex p-power roots
of unity. Haar measure is real, normalized by vol(Z_p)=1. These integrals
are distinct from the Q_p-valued measures in transforms and measures.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from itertools import islice

from ._cyclotomic import CyclotomicValue, _size
from .core import DomainError, Padic, PadicField, PrecisionError, exact, is_prime, p_power, valuation


def _prime(field=None, value=None):
    if field is None and isinstance(value,Padic):
        return value.p
    p = field.p if isinstance(field,PadicField) else field
    if p is None:
        raise TypeError("exact arguments require a prime or PadicField")
    if not is_prime(p):
        raise ValueError("p must be prime")
    if isinstance(value,Padic) and value.p != p:
        raise TypeError("p-adic argument has a different prime")
    return p


def additive_phase(value, field=None):
    """The exact phase {x}_p in [0,1) modulo Z_p for chi(x)=exp(2*pi*i*{x}_p).

    Accept exact rational input with a prime/field, or infer the prime from
    a Padic ball. A ball needs absolute precision >=0: otherwise its phase
    is not determined. This is the p-adic fractional part, not the real one.
    """
    p = _prime(field,value)
    if isinstance(value,Padic):
        if value.precision < 0:
            raise PrecisionError("additive character requires absolute precision >=0")
        value = value.representative()
    q = exact(value)
    v = valuation(q,p)
    if v >= 0:
        return Fraction(0)
    denominator = p**(-v)
    integral = q*denominator
    residue = integral.numerator*pow(integral.denominator,-1,denominator) % denominator
    return Fraction(residue,denominator)


def additive_character(value, field=None):
    """Exact complex root chi_p(x), with kernel Z_p and positive phase.

    The root order is limited to 4096 by the exact coefficient backend.
    additive_phase itself has no such conductor limit. Use .complex(dps)
    on the returned CyclotomicValue for an archimedean approximation.
    """
    p = _prime(field,value)
    phase = additive_phase(value,p)
    denominator,level = phase.denominator,0
    while denominator > 1:
        denominator //= p
        level += 1
    return CyclotomicValue.root(p,level,phase.numerator)


def haar_volume(field, exponent=0):
    """Exact real Haar volume of any a+p**exponent*Z_p: p**(-exponent)."""
    p = _prime(field)
    if type(exponent) is not int:
        raise TypeError("ball exponent must be an integer")
    return p_power(p,-exponent)


def _values(values, p, max_size):
    _size(p,0,max_size)
    values = tuple(islice(iter(values),max_size+1))
    if not values:
        raise ValueError("at least one function value is required")
    if len(values) > max_size:
        raise PrecisionError(f"number of values exceeds max_size={max_size}")
    result = []
    for value in values:
        if isinstance(value,CyclotomicValue):
            if value.p != p:
                if value.level:
                    raise TypeError("function coefficients have a different cyclotomic prime")
                value = CyclotomicValue(p,[value.as_fraction()])
        else:
            value = CyclotomicValue(p,[exact(value)])
        result.append(value)
    return tuple(result)


def finite_fourier(values, field, *, inverse=False, max_size=256):
    """Normalized forward DFT on Z/(p**n), with an unnormalized inverse.

    Forward coefficient k is sum_j(f_j*exp(-2*pi*i*j*k/N))/N, N=p**n.
    Input and output coefficients are exact CyclotomicValue objects. The
    naive native algorithm is intended for small finite models. max_size
    caps N; the coefficient backend additionally caps root order at 4096.
    """
    p = _prime(field)
    if type(inverse) is not bool:
        raise TypeError("inverse must be a bool")
    values = _values(values,p,max_size)
    size,n = len(values),0
    remaining = size
    while remaining > 1 and remaining % p == 0:
        remaining //= p
        n += 1
    if remaining != 1:
        raise ValueError("number of values must be a power of p, including p**0")
    level = max(n,*(value.level for value in values))
    order = _size(p,level,4096)
    step = order//size
    embedded = [value._terms(level) for value in values]
    sign,scale = (1,Fraction(1)) if inverse else (-1,Fraction(1,size))
    result = []
    for k in range(size):
        terms = {}
        for j,coefficients in enumerate(embedded):
            shift = sign*j*k*step
            for i,c in coefficients.items():
                index = (i+shift) % order
                terms[index] = terms.get(index,Fraction(0))+c*scale
        result.append(CyclotomicValue._from_terms(p,level,terms))
    return tuple(result)


@dataclass(frozen=True, init=False)
class PadicTestFunction:
    """A compactly supported, locally constant function with exact values.

    Values index cosets p**a*j + p**b*Z_p in p**a*Z_p, with
    a=support_exponent <= b=constancy_exponent and N=p**(b-a) values.
    The function vanishes outside that container. Values are rational or
    CyclotomicValue coefficients, with a fixed complex root-of-unity place.
    Field precision is not a sampling precision; this is complete finite
    data for a function on Q_p. max_size bounds its number of cosets.
    """
    field: PadicField
    values: tuple[CyclotomicValue, ...]
    support_exponent: int
    constancy_exponent: int

    def __init__(self, field, values, *, support_exponent=0,
                 constancy_exponent=1, max_size=4096):
        p = _prime(field)
        a,b = support_exponent,constancy_exponent
        if type(a) is not int or type(b) is not int:
            raise TypeError("support and constancy exponents must be integers")
        if b < a:
            raise ValueError("constancy exponent must be >= support exponent")
        size = _size(p,b-a,max_size)
        values = _values(values,p,max_size)
        if len(values) != size:
            raise ValueError(f"this lattice needs {size} values")
        for key,value in (("field",field if isinstance(field,PadicField) else PadicField(p)),
                          ("values",values),("support_exponent",a),("constancy_exponent",b)):
            object.__setattr__(self,key,value)

    @classmethod
    def ball(cls, field, center=0, exponent=0, *, max_size=4096):
        """The exact indicator of center+p**exponent*Z_p.

        A noncentral ball is represented inside a sufficiently large origin
        container. Input centre precision must determine the ball.
        """
        p = _prime(field,center)
        if type(exponent) is not int:
            raise TypeError("ball exponent must be an integer")
        if isinstance(center,Padic):
            if center.precision < exponent:
                raise PrecisionError("centre digits do not determine the requested ball")
            center = center.representative()
        center = exact(center)
        a = min(exponent,valuation(center,p))
        size = _size(p,exponent-a,max_size)
        index = cls._residue(center,p,a,exponent)
        return cls(field,[int(j == index) for j in range(size)],
                   support_exponent=a,constancy_exponent=exponent,max_size=max_size)

    @staticmethod
    def _residue(point, p, a, b):
        scaled = point/p_power(p,a)
        modulus = p**(b-a)
        return scaled.numerator*pow(scaled.denominator,-1,modulus) % modulus

    def _index(self, point):
        p,a,b = self.field.p,self.support_exponent,self.constancy_exponent
        _prime(p,point)
        if isinstance(point,Padic):
            if point.unit and point.v < a:
                return None
            if point.precision < b:
                raise PrecisionError("point digits do not determine the function coset")
            point = point.representative()
        point = exact(point)
        if valuation(point,p) < a:
            return None
        return self._residue(point,p,a,b)

    def __call__(self, point):
        """Exact evaluation; a p-adic ball must determine the required coset."""
        index = self._index(point)
        return CyclotomicValue(self.field.p) if index is None else self.values[index]

    def _replace_values(self, values):
        return type(self)(self.field,values,support_exponent=self.support_exponent,
                          constancy_exponent=self.constancy_exponent,
                          max_size=max(4096,len(self.values)))

    def _check(self, other):
        if not isinstance(other,PadicTestFunction):
            raise TypeError("expected a PadicTestFunction")
        if (self.field.p,self.support_exponent,self.constancy_exponent) != (
                other.field.p,other.support_exponent,other.constancy_exponent):
            raise DomainError("function operation requires the same prime and lattice")

    def integral(self):
        """Exact Haar integral p**(-b)*sum(values)."""
        return sum(self.values,CyclotomicValue(self.field.p))*haar_volume(
            self.field,self.constancy_exponent)

    def inner_product(self, other):
        """Exact integral self(x)*conj(other(x)) dx on a common lattice."""
        self._check(other)
        return sum((a*b.conjugate() for a,b in zip(self.values,other.values)),
                   CyclotomicValue(self.field.p))*haar_volume(self.field,self.constancy_exponent)

    def norm_squared(self):
        """Exact Haar L2 norm squared, with the chosen complex conjugation."""
        return self.inner_product(self)

    def __mul__(self, other):
        if isinstance(other,PadicTestFunction):
            self._check(other)
            return self._replace_values(a*b for a,b in zip(self.values,other.values))
        scalar = CyclotomicValue(self.field.p)._coerce(other)
        return self._replace_values(value*scalar for value in self.values)

    __rmul__ = __mul__

    def __add__(self, other):
        self._check(other)
        return self._replace_values(a+b for a,b in zip(self.values,other.values))

    def reflect(self):
        """The exact function x -> self(-x)."""
        size = len(self.values)
        return self._replace_values(self.values[-j % size] for j in range(size))

    def translate(self, shift):
        """x -> self(x-shift), for a shift inside the support container."""
        index = self._index(shift)
        if index is None:
            raise DomainError("translation must stay inside the support container")
        size = len(self.values)
        return self._replace_values(self.values[(j-index) % size] for j in range(size))

    def fourier(self, *, inverse=False, max_size=256):
        """Exact Q_p Fourier transform with kernel chi(-x*xi).

        The dual lattice has support exponent -b and constancy exponent -a.
        Haar weights are p**(-b); inverse=True uses the positive kernel.
        The finite data describes the entire transform, without an analytic
        cutoff. A default limit of 256 cosets bounds the native DFT work.
        """
        p,a,b = self.field.p,self.support_exponent,self.constancy_exponent
        values = finite_fourier(self.values,p,inverse=inverse,max_size=max_size)
        factor = p_power(p,-b if inverse else -a)
        return type(self)(self.field,[value*factor for value in values],
                          support_exponent=-b,constancy_exponent=-a,
                          max_size=max_size)

    def convolve(self, other, *, max_size=256):
        """Exact additive Haar convolution on a common lattice.

        Result k is p**(-b)*sum_j(self[j]*other[k-j]). The convolution
        theorem uses these real Haar weights, rather than a Q_p measure.
        """
        self._check(other)
        _size(self.field.p,0,max_size)
        size = len(self.values)
        if size > max_size:
            raise PrecisionError(f"convolution needs {size} cosets; max_size={max_size}")
        factor = haar_volume(self.field,self.constancy_exponent)
        values = [sum((a*other.values[(k-j) % size] for j,a in enumerate(self.values)),
                      CyclotomicValue(self.field.p))*factor for k in range(size)]
        return self._replace_values(values)
