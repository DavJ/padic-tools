"""Unramified quadratic fields, exact residues and shared ideal-ball arithmetic.

The integral power basis is certified by irreducibility modulo p. Every
local element is one ideal ball a + b*theta + p**N O_L, not a formal quotient
with an unverified valuation. No optional computer algebra system is used.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from .core import DomainError, Padic, PadicField, PrecisionError, exact, is_prime, p_power, valuation
from .symbolic import Polynomial


def _mod_p(value, p):
    q = exact(value)
    if q.denominator % p == 0:
        raise DomainError("reduction modulo p requires a p-integral rational")
    return q.numerator * pow(q.denominator, -1, p) % p


def _defining_polynomial(p, polynomial):
    if not is_prime(p):
        raise ValueError("p must be prime")
    f = polynomial if isinstance(polynomial, Polynomial) else Polynomial(polynomial)
    if f.degree != 2 or f[2] != 1:
        raise DomainError("the defining polynomial must be monic of degree two")
    c, b = _mod_p(f[0], p), _mod_p(f[1], p)
    if p == 2:
        irreducible = c != 0 and (1 + b + c) % p != 0
    else:
        discriminant = (b*b - 4*c) % p
        irreducible = pow(discriminant, (p-1)//2, p) == p-1
    if not irreducible:
        raise DomainError("the defining polynomial must be irreducible modulo p; "
                          "split and ramified presentations are unsupported")
    return f


def _name(name):
    if not isinstance(name, str) or not name.isidentifier():
        raise ValueError("the generator name must be an identifier")
    return name


def _multiply(x, y, b, c):
    a0, a1 = x
    d0, d1 = y
    return a0*d0 - c*a1*d1, a0*d1 + a1*d0 - b*a1*d1


@dataclass(frozen=True, init=False, eq=False)
class QuadraticResidueField:
    """Exact F_p[t]/(f), with a checked irreducible monic quadratic f.

    Coefficients of ``polynomial`` are in ascending degree order. Normally
    obtain this object from ``UnramifiedQuadraticField.residue_field()``.
    """
    p: int
    modulus: tuple[int, int, int]
    name: str

    def __init__(self, p, polynomial, name="a"):
        f = _defining_polynomial(p, polynomial)
        object.__setattr__(self, "p", p)
        object.__setattr__(self, "modulus", (_mod_p(f[0], p), _mod_p(f[1], p), 1))
        object.__setattr__(self, "name", _name(name))

    @property
    def order(self):
        """Cardinality p**2."""
        return self.p**2

    @property
    def generator(self):
        """The exact residue class of t."""
        return self(0, 1)

    def __call__(self, value=0, b=None):
        if isinstance(value, QuadraticResidueElement):
            if b is not None or value.field != self:
                raise TypeError("residue fields must have the same prime and modulus")
            return QuadraticResidueElement(self, value.a, value.b)
        if b is None:
            if isinstance(value, (tuple, list)):
                if len(value) != 2:
                    raise ValueError("two coefficients are required")
                value, b = value
            else:
                b = 0
        return QuadraticResidueElement(self, _mod_p(value, self.p), _mod_p(b, self.p))

    def __eq__(self, other):
        if not isinstance(other, QuadraticResidueField):
            return NotImplemented
        return (self.p, self.modulus) == (other.p, other.modulus)

    def __hash__(self):
        return hash((self.p, self.modulus))

    def __repr__(self):
        return f"F_{self.p}^2({self.name}; modulus={self.modulus})"


@dataclass(frozen=True)
class QuadraticResidueElement:
    """An exact pair of residues in a QuadraticResidueField."""
    field: QuadraticResidueField
    a: int
    b: int

    def __post_init__(self):
        if not isinstance(self.field, QuadraticResidueField):
            raise TypeError("an exact quadratic residue field is required")
        if type(self.a) is not int or type(self.b) is not int:
            raise TypeError("residue coefficients must be integers")
        object.__setattr__(self, "a", self.a % self.field.p)
        object.__setattr__(self, "b", self.b % self.field.p)

    @property
    def coefficients(self):
        """Canonical coefficients (a,b), each in range(p)."""
        return self.a, self.b

    def __add__(self, other):
        try:
            other = self.field(other)
        except TypeError:
            return NotImplemented
        return self.field(self.a + other.a, self.b + other.b)

    __radd__ = __add__

    def __neg__(self):
        return self.field(-self.a, -self.b)

    def __sub__(self, other):
        try:
            return self + (-self.field(other))
        except TypeError:
            return NotImplemented

    def __rsub__(self, other):
        return -self + other

    def __mul__(self, other):
        try:
            other = self.field(other)
        except TypeError:
            return NotImplemented
        c, b, _ = self.field.modulus
        return self.field(_multiply(self.coefficients, other.coefficients, b, c))

    __rmul__ = __mul__

    def conjugate(self):
        """The nontrivial automorphism, also the p-power Frobenius."""
        return self.field(self.a - self.field.modulus[1]*self.b, -self.b)

    def frobenius(self):
        """The exact p-power automorphism."""
        return self.conjugate()

    def trace(self):
        """Trace to F_p, returned as an integer residue."""
        return (2*self.a - self.field.modulus[1]*self.b) % self.field.p

    def norm(self):
        """Norm to F_p, returned as an integer residue."""
        c, b, _ = self.field.modulus
        return (self.a*self.a - b*self.a*self.b + c*self.b*self.b) % self.field.p

    def inverse(self):
        """Exact inverse; zero raises ZeroDivisionError."""
        if not self:
            raise ZeroDivisionError("zero has no inverse in the residue field")
        return self.conjugate() * pow(self.norm(), -1, self.field.p)

    def __truediv__(self, other):
        try:
            return self * self.field(other).inverse()
        except TypeError:
            return NotImplemented

    def __rtruediv__(self, other):
        return self.inverse() * other

    def __pow__(self, exponent):
        if type(exponent) is not int:
            raise TypeError("the residue-field exponent must be an integer")
        if exponent < 0:
            return self.inverse()**(-exponent)
        result, base = self.field(1), self
        while exponent:
            if exponent & 1:
                result *= base
            exponent >>= 1
            if exponent:
                base *= base
        return result

    def __bool__(self):
        return self.a != 0 or self.b != 0

    def __repr__(self):
        return f"{self.a} + {self.b}*{self.field.name} (mod {self.field.p})"


@dataclass(frozen=True)
class _QuadraticField:
    """Shared power-basis hooks; public subclasses certify their basis."""
    base: PadicField
    polynomial: Polynomial
    name: str = "a"

    def __post_init__(self):
        if not isinstance(self.base, PadicField):
            raise TypeError("the base must be a Qp field")
        object.__setattr__(self, "polynomial", _defining_polynomial(self.p, self.polynomial))
        object.__setattr__(self, "name", _name(self.name))

    @property
    def p(self):
        """Residue characteristic and chosen uniformizer."""
        return self.base.p

    @property
    def precision(self):
        """Default relative precision in uniformizer digits for exact inputs."""
        return self.ramification_index*self.base.precision

    @property
    def degree(self):
        """Extension degree over the base field: two."""
        return 2

    @property
    def ramification_index(self):
        """Ramification index: one."""
        return 1

    @property
    def residue_degree(self):
        """Degree of the residue field over F_p: two."""
        return 2

    @property
    def generator(self):
        """The chosen root theta, as a ball at default relative precision."""
        return self(0, 1)

    def uniformizer(self):
        """The element p, with valuation one."""
        return self(self.p)

    def residue_field(self):
        """The exact F_(p**2) in the matching reduced power basis."""
        return QuadraticResidueField(self.p, self.polynomial, self.name + "0")

    def lift_residue(self, value, *, absprec=None):
        """Choose integer coefficient representatives; this is not a ring map."""
        r = self.residue_field()(value)
        n = max(1, self.precision) if absprec is None else absprec
        if type(n) is not int or n < 1:
            raise PrecisionError("a residue lift needs absolute precision at least one")
        return self(*r.coefficients, absprec=n)

    def _weights(self):
        return 0, 0

    def _valuation(self, coordinates):
        return min(self.ramification_index*valuation(q, self.p)+w
                   for q, w in zip(coordinates, self._weights()))

    def _coordinate_precisions(self, n):
        if type(n) is not int:
            raise TypeError("absolute precision must be an integer")
        e = self.ramification_index
        return tuple((n-w+e-1)//e for w in self._weights())

    def _trace_precision(self, n):
        return n

    def _norm_precision(self, n, v):
        return n+v

    def _residue(self, element):
        return self.residue_field()(element.a.residue(1), element.b.residue(1))

    def _new_element(self, coordinates, n):
        return UnramifiedQuadraticElement._from_pair(self, coordinates, n)

    def _uniformizer_name(self):
        return str(self.p)

    def _key(self):
        return type(self), self.p, self.polynomial

    def _data(self, value, b=None):
        if isinstance(value, _QuadraticElement):
            if b is not None or value.field._key() != self._key():
                raise TypeError("extension presentations must have the same prime and polynomial")
            return value.representative(), value.precision
        if b is None:
            if isinstance(value, (tuple, list)):
                if len(value) != 2:
                    raise ValueError("two coefficients are required")
                value, b = value
            else:
                b = 0
        coordinates, precisions = [], []
        for coefficient, weight in zip((value, b), self._weights()):
            if isinstance(coefficient, Padic):
                if coefficient.p != self.p:
                    raise TypeError("the coefficient prime differs from the extension prime")
                coordinates.append(coefficient.representative())
                precisions.append(self.ramification_index*coefficient.precision+weight)
            else:
                coordinates.append(exact(coefficient))
        return tuple(coordinates), min(precisions) if precisions else None

    def __call__(self, value=0, b=None, *, absprec=None):
        coordinates, known = self._data(value, b)
        if absprec is None:
            v = self._valuation(coordinates)
            n = known if known is not None else self.precision + (v if any(coordinates) else 0)
        else:
            n = absprec
        if type(n) is not int:
            raise TypeError("absolute precision must be an integer")
        if known is not None and n > known:
            raise PrecisionError("cannot manufacture digits of an existing ball")
        return self._new_element(coordinates, n)

    def O(self, exponent):
        """The uniformizer ideal ball, including negative exponents."""
        return self._new_element((Fraction(0), Fraction(0)), exponent)

    def __repr__(self):
        return f"{type(self).__name__}({self.base!r}, {self.polynomial.coefficients!r}, name={self.name!r})"


@dataclass(frozen=True, repr=False)
class UnramifiedQuadraticField(_QuadraticField):
    """Q_p[t]/(f), for monic p-integral quadratic f irreducible modulo p.

    ``base`` is a Qp field. ``polynomial`` is a Polynomial or ascending
    exact rational coefficients. ``L(a,b)`` creates a+b*theta with the base
    field's default relative precision; ``absprec=N`` specifies an ideal
    ball modulo p**N. Existing balls can only lose digits on conversion.
    Different polynomial presentations require an explicit isomorphism.
    """


@dataclass(frozen=True, init=False, eq=False)
class _QuadraticElement:
    """Shared arithmetic, with valuation and precision supplied by the field."""
    field: _QuadraticField
    a: Padic
    b: Padic
    _precision: int

    def __new__(cls, *args, **kwargs):
        raise TypeError("construct extension balls through a quadratic local field")

    @classmethod
    def _from_pair(cls, field, coordinates, n):
        na, nb = field._coordinate_precisions(n)
        obj = object.__new__(cls)
        object.__setattr__(obj, "field", field)
        object.__setattr__(obj, "a", field.base(coordinates[0], absprec=na))
        object.__setattr__(obj, "b", field.base(coordinates[1], absprec=nb))
        object.__setattr__(obj, "_precision", n)
        return obj

    @property
    def p(self):
        """Residue characteristic."""
        return self.field.p

    @property
    def precision(self):
        """Absolute precision N in powers of the field's uniformizer."""
        return self._precision

    @property
    def v(self):
        """Valuation, or its lower bound N for a zero ball."""
        return min(self.field.ramification_index*q.v+w
                   for q, w in zip(self.coefficients, self.field._weights()))

    @property
    def relative_precision(self):
        """Number of known significant digits, zero for a zero ball."""
        return self.precision - self.v

    @property
    def coefficients(self):
        """Qp coordinate balls with basis-dependent absolute precisions."""
        return self.a, self.b

    def representative(self):
        """One rational coordinate pair, not an embedding into C."""
        return self.a.representative(), self.b.representative()

    def is_zero(self):
        """Whether the ideal ball contains zero, not exact equality to zero."""
        return self.a.is_zero() and self.b.is_zero()

    def valuation(self):
        """Known valuation with v(uniformizer)=1; zero balls raise."""
        if self.is_zero():
            raise PrecisionError(f"only valuation >= {self.precision} is known")
        return self.v

    def base_valuation(self):
        """Valuation normalized by v(p)=1, as an exact rational."""
        return Fraction(self.valuation(), self.field.ramification_index)

    def absolute_value_squared(self):
        """Exact square of the absolute value extending the base p-adic norm."""
        return p_power(self.p, -2*self.valuation()//self.field.ramification_index)

    def contains(self, value, b=None):
        """Test membership of an exact rational coordinate pair or scalar."""
        coordinates, known = self.field._data(value, b)
        if known is not None:
            raise TypeError("membership requires exact coordinates; use is_congruent for balls")
        return self.field._valuation(tuple(x-y for x, y in
                                          zip(self.representative(), coordinates))) >= self.precision

    def is_congruent(self, other, precision=None):
        """Compare known coordinates modulo p**precision; never invent digits."""
        coordinates, known = self.field._data(other)
        known = self.precision if known is None else min(self.precision, known)
        n = known if precision is None else precision
        if type(n) is not int:
            raise TypeError("precision must be an integer")
        if n > known:
            raise PrecisionError("requested congruence exceeds known precision")
        return self.field._valuation(tuple(x-y for x, y in
                                          zip(self.representative(), coordinates))) >= n

    def truncate(self, precision):
        """Forget digits of both coordinates, never increase precision."""
        return self.field(self, absprec=precision)

    def residue(self):
        """Reduce an integral ball to its field's exact residue field."""
        if self.v < 0:
            raise DomainError("a nonintegral extension element has no residue")
        if self.precision < 1:
            raise PrecisionError("reduction needs at least one known uniformizer digit")
        return self.field._residue(self)

    def __add__(self, other):
        try:
            coordinates, known = self.field._data(other)
        except TypeError:
            return NotImplemented
        n = self.precision if known is None else min(self.precision, known)
        return self._from_pair(self.field, tuple(x+y for x, y in zip(self.representative(), coordinates)), n)

    __radd__ = __add__

    def __neg__(self):
        return self._from_pair(self.field, tuple(-q for q in self.representative()), self.precision)

    def __sub__(self, other):
        try:
            coordinates, known = self.field._data(other)
        except TypeError:
            return NotImplemented
        n = self.precision if known is None else min(self.precision, known)
        return self._from_pair(self.field, tuple(x-y for x, y in zip(self.representative(), coordinates)), n)

    def __rsub__(self, other):
        return -self + other

    def __mul__(self, other):
        try:
            coordinates, known = self.field._data(other)
        except TypeError:
            return NotImplemented
        v = self.field._valuation(coordinates)
        if known is None:
            n = self.precision + (v if any(coordinates) else 0)
        else:
            v = min(v, known)
            n = min(self.precision + v, known + self.v)
        f = self.field.polynomial
        return self._from_pair(self.field, _multiply(self.representative(), coordinates, f[1], f[0]), n)

    __rmul__ = __mul__

    def conjugate(self):
        """The nontrivial field automorphism; it preserves absolute precision."""
        a, b = self.representative()
        return self._from_pair(self.field, (a-self.field.polynomial[1]*b, -b), self.precision)

    def trace(self):
        """Algebraic trace to Qp, with the certified trace-ideal precision."""
        a, b = self.representative()
        return self.field.base(2*a-self.field.polynomial[1]*b,
                               absprec=self.field._trace_precision(self.precision))

    def norm(self):
        """Algebraic norm to Qp, with the field's norm-error bound."""
        a, b = self.representative()
        f = self.field.polynomial
        return self.field.base(a*a-f[1]*a*b+f[0]*b*b,
                               absprec=self.field._norm_precision(self.precision, self.v))

    def inverse(self):
        """Invert a ball excluding zero, with absolute precision N-2*v."""
        if self.is_zero():
            raise PrecisionError("cannot invert an extension ball containing zero")
        a, b = self.representative()
        f = self.field.polynomial
        norm = a*a-f[1]*a*b+f[0]*b*b
        return self._from_pair(self.field, ((a-f[1]*b)/norm, -b/norm), self.precision-2*self.v)

    def __truediv__(self, other):
        try:
            coordinates, known = self.field._data(other)
        except TypeError:
            return NotImplemented
        if known is not None:
            return self * self.field(other).inverse()
        a, b = coordinates
        f = self.field.polynomial
        norm = a*a-f[1]*a*b+f[0]*b*b
        if not norm:
            raise ZeroDivisionError("division by exact zero")
        return self * ((a-f[1]*b)/norm, -b/norm)

    def __rtruediv__(self, other):
        return self.inverse() * other

    def __pow__(self, exponent):
        if type(exponent) is not int:
            raise TypeError("the extension-field exponent must be an integer")
        if exponent < 0:
            return self.inverse()**(-exponent)
        if exponent == 0:
            return self.field(1)
        result, base = None, self
        while exponent:
            if exponent & 1:
                result = base if result is None else result*base
            exponent >>= 1
            if exponent:
                base *= base
        return result

    def __bool__(self):
        raise TypeError("extension-ball truth is ambiguous; use is_zero() or contains()")

    def __eq__(self, other):
        if not isinstance(other, _QuadraticElement):
            return NotImplemented
        return (self.field._key(), self.representative(), self.precision) == (
            other.field._key(), other.representative(), other.precision)

    def __hash__(self):
        return hash((self.field._key(), self.representative(), self.precision))

    def __repr__(self):
        a, b = self.representative()
        return f"{a} + {b}*{self.field.name} + O({self.field._uniformizer_name()}^{self.precision})"


@dataclass(frozen=True, init=False, eq=False, repr=False)
class UnramifiedQuadraticElement(_QuadraticElement):
    """An immutable ideal ball a+b*theta + O(p**N) in a checked field.

    ``a`` and ``b`` are Qp coordinate balls with the same absolute precision.
    ``norm()`` is the algebraic norm to Qp; ``absolute_value()`` is p**(-v).
    Equality compares identical balls. ``is_congruent`` compares known digits.
    """
    field: UnramifiedQuadraticField

    def absolute_value(self):
        """The absolute value p**(-v), for a ball excluding zero."""
        return p_power(self.p, -self.valuation())

    def frobenius(self):
        """The lift of residue p-power Frobenius; equals conjugation in degree two."""
        return self.conjugate()
