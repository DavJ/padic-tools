"""Quadratic Eisenstein fields, with uniformizer precision and exact residues.

Q2(i) is presented in the integral basis 1, pi for pi=1+i, whose polynomial
is X**2-2*X+2. The Gaussian basis 1,i is not a weighted valuation basis.
"""
from __future__ import annotations

from dataclasses import dataclass

from .core import DomainError, PadicField, PrecisionError, Qp, is_prime, valuation
from .extensions import _QuadraticElement, _QuadraticField, _mod_p, _name
from .symbolic import Polynomial


@dataclass(frozen=True)
class PrimeResidueField:
    """The exact prime field F_p, used as an Eisenstein field's residue field."""
    p: int

    def __post_init__(self):
        if not is_prime(self.p):
            raise ValueError("p must be prime")

    @property
    def order(self):
        """Cardinality p."""
        return self.p

    def __call__(self, value=0):
        if isinstance(value, PrimeResidueElement):
            if value.field != self:
                raise TypeError("residue fields must have the same prime")
            return value
        return PrimeResidueElement(self, _mod_p(value, self.p))

    def __repr__(self):
        return f"F_{self.p}"


@dataclass(frozen=True)
class PrimeResidueElement:
    """An exact integer residue in a PrimeResidueField."""
    field: PrimeResidueField
    value: int

    def __post_init__(self):
        if not isinstance(self.field, PrimeResidueField):
            raise TypeError("an exact prime residue field is required")
        if type(self.value) is not int:
            raise TypeError("the residue coefficient must be an integer")
        object.__setattr__(self, "value", self.value % self.field.p)

    def __add__(self, other):
        try:
            other = self.field(other)
        except TypeError:
            return NotImplemented
        return self.field(self.value+other.value)

    __radd__ = __add__

    def __neg__(self):
        return self.field(-self.value)

    def __sub__(self, other):
        try:
            return self + (-self.field(other))
        except TypeError:
            return NotImplemented

    def __rsub__(self, other):
        return -self+other

    def __mul__(self, other):
        try:
            other = self.field(other)
        except TypeError:
            return NotImplemented
        return self.field(self.value*other.value)

    __rmul__ = __mul__

    def inverse(self):
        """Exact inverse; zero raises ZeroDivisionError."""
        if not self.value:
            raise ZeroDivisionError("zero has no inverse in the residue field")
        return self.field(pow(self.value, -1, self.field.p))

    def __truediv__(self, other):
        try:
            return self*self.field(other).inverse()
        except TypeError:
            return NotImplemented

    def __rtruediv__(self, other):
        return self.inverse()*other

    def __pow__(self, exponent):
        if type(exponent) is not int:
            raise TypeError("the residue-field exponent must be an integer")
        if exponent < 0:
            return self.inverse()**(-exponent)
        return self.field(pow(self.value, exponent, self.field.p))

    def frobenius(self):
        """Frobenius on the prime field is the identity."""
        return self

    def trace(self):
        """Trace of F_p over itself, as an integer residue."""
        return self.value

    def norm(self):
        """Norm of F_p over itself, as an integer residue."""
        return self.value

    def __bool__(self):
        return self.value != 0

    def __repr__(self):
        return f"{self.value} (mod {self.field.p})"


@dataclass(frozen=True, repr=False)
class EisensteinQuadraticField(_QuadraticField):
    """Q_p[pi] for a certified monic Eisenstein quadratic X**2+B*X+C.

    Require v_p(B)>=1 and v_p(C)=1. The generator pi is a uniformizer,
    with v_pi(p)=2 and residue field F_p. ``L(a,b,absprec=N)`` represents
    a+b*pi modulo pi**N; coordinate precisions are ceil(N/2), floor(N/2).
    The default relative precision is twice ``base.precision``. Explicit
    precisions and valuations are always in uniformizer units.
    """
    name: str = "pi"

    def __post_init__(self):
        if not isinstance(self.base, PadicField):
            raise TypeError("the base must be a Qp field")
        f = self.polynomial
        f = f if isinstance(f, Polynomial) else Polynomial(f)
        if f.degree != 2 or f[2] != 1:
            raise DomainError("the defining polynomial must be monic of degree two")
        if valuation(f[1], self.p) < 1 or valuation(f[0], self.p) != 1:
            raise DomainError("the polynomial must be Eisenstein: v_p(B)>=1 and v_p(C)=1")
        object.__setattr__(self, "polynomial", f)
        object.__setattr__(self, "name", _name(self.name))

    @property
    def ramification_index(self):
        """Ramification index: two."""
        return 2

    @property
    def residue_degree(self):
        """Residue degree: one."""
        return 1

    @property
    def different_exponent(self):
        """v_pi(f'(pi)), the exponent of the different ideal."""
        return min(2*valuation(self.polynomial[1], self.p),
                   2*valuation(2, self.p)+1)

    def _weights(self):
        return 0, 1

    def uniformizer(self):
        """The generator pi, of valuation one."""
        return self.generator

    def residue_field(self):
        """The exact prime field F_p."""
        return PrimeResidueField(self.p)

    def lift_residue(self, value, *, absprec=None):
        """Choose the integer coefficient lift; this is not a ring map."""
        r = self.residue_field()(value)
        n = self.precision if absprec is None else absprec
        if type(n) is not int or n < 1:
            raise PrecisionError("a residue lift needs absolute precision at least one")
        return self(r.value, absprec=n)

    def _trace_precision(self, n):
        return (n+self.different_exponent)//2

    def _norm_precision(self, n, v):
        return min(self._trace_precision(n+v), n)

    def _residue(self, element):
        return self.residue_field()(element.a.residue(1))

    def _new_element(self, coordinates, n):
        return EisensteinQuadraticElement._from_pair(self, coordinates, n)

    def _uniformizer_name(self):
        return self.name


@dataclass(frozen=True, init=False, eq=False, repr=False)
class EisensteinQuadraticElement(_QuadraticElement):
    """A ball a+b*pi + O(pi**N) in a quadratic Eisenstein field.

    Its Qp coordinates have absolute precisions ceil(N/2), floor(N/2).
    ``valuation()`` and ``precision`` use pi units; ``base_valuation()``
    uses v(p)=1. ``norm()`` and ``trace()`` return Qp balls with p precision.
    ``absolute_value_squared()`` gives an exact rational magnitude, including
    when the absolute value itself is an irrational square root of p.
    """
    field: EisensteinQuadraticField


def Q2i(precision=20):
    """Q2(i), presented with uniformizer pi=1+i and X**2-2*X+2.

    ``precision`` counts base-2 significant digits; the returned field's
    default relative precision is twice that many pi digits. Obtain i as
    ``L.generator-1``. Coordinates supplied to L use the basis (1,pi).
    """
    return EisensteinQuadraticField(Qp(2, precision), [2, -2, 1])
