"""Finite precision balls in Q_p, with explicit absolute error bounds.

Every Padic denotes a ball a + p**N Z_p, never a hidden exact rational.
Exact rational scalars in arithmetic are handled without rounding them first.
"""
from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal
from fractions import Fraction
from functools import lru_cache
from math import inf


class PrecisionError(ArithmeticError):
    """The known digits do not justify the requested operation or conclusion."""


class DomainError(ValueError):
    """An argument is outside an explicitly supported mathematical domain."""


def exact(value) -> Fraction:
    """Convert int, Fraction, Decimal or decimal/rational text; reject floats.

    Use Fraction.from_float explicitly if binary floating point was intended.
    """
    if isinstance(value, (float, complex, bool)):
        raise TypeError("use int, Fraction, Decimal or text, not float/complex/bool")
    if isinstance(value, (int, Fraction, Decimal, str)):
        return Fraction(value)
    raise TypeError(f"unsupported exact scalar: {type(value).__name__}")


@lru_cache(maxsize=256)
def is_prime(p: int) -> bool:
    """Deterministic primality test for integers below 2**64.

    Uses the deterministic seven-base Miller--Rabin test in this range.
    Larger integers require a different certified primality backend.
    """
    if type(p) is not int:
        raise TypeError("p must be an integer")
    if p >= 2**64:
        raise ValueError("native backend supports certified primes below 2**64")
    if p < 2:
        return False
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if p % q == 0:
            return p == q
    d, r = p - 1, 0
    while d % 2 == 0:
        d //= 2
        r += 1
    for a in (2, 325, 9375, 28178, 450775, 9780504, 1795265022):
        if a % p == 0:
            continue
        y = pow(a, d, p)
        if y in (1, p - 1):
            continue
        for _ in range(r - 1):
            y = y * y % p
            if y == p - 1:
                break
        else:
            return False
    return True


def _vp_integer(n: int, p: int) -> int:
    n = abs(n)
    if not n:
        raise ValueError("zero has infinite exact valuation")
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v


def valuation(value, p: int):
    """Exact p-adic valuation of a rational (infinity for exact zero)."""
    if not is_prime(p):
        raise ValueError("p must be prime")
    x = exact(value)
    if not x:
        return inf
    return _vp_integer(x.numerator, p) - _vp_integer(x.denominator, p)


def p_power(p: int, exponent: int) -> Fraction:
    """An exact rational power of p, allowing negative exponents."""
    return Fraction(p**exponent) if exponent >= 0 else Fraction(1, p**(-exponent))


@dataclass(frozen=True)
class PadicField:
    """Q_p with a default number of significant digits for exact inputs.

    ``Qp(5, 20)(25)`` has absolute precision 22. Explicit ``absprec``
    overrides that default for exact input, but cannot extend an existing ball.
    """
    p: int
    precision: int = 20

    def __post_init__(self):
        if not is_prime(self.p):
            raise ValueError("p must be prime")
        if type(self.precision) is not int or self.precision < 1:
            raise ValueError("precision must be a positive integer")

    def __call__(self, value=0, *, absprec=None) -> Padic:
        if isinstance(value, Padic):
            if self.p != value.p:
                raise TypeError("cannot coerce between different p-adic primes")
            n = value.precision if absprec is None else absprec
            if n > value.precision:
                raise PrecisionError("cannot manufacture digits of an existing ball")
            return Padic._from_fraction(self, value.representative(), n)
        x = exact(value)
        n = (self.precision + (valuation(x, self.p) if x else 0)
             if absprec is None else absprec)
        return Padic._from_fraction(self, x, n)

    def O(self, exponent: int) -> Padic:
        """The zero ball O(p**exponent), with unknown higher digits."""
        return Padic._from_fraction(self, Fraction(0), exponent)

    def __repr__(self):
        return f"Qp({self.p}, precision={self.precision})"


Qp = PadicField


@dataclass(frozen=True, eq=False, init=False)
class Padic:
    """An immutable ball ``unit * p**v + O(p**precision)``.

    Construct through a PadicField. Equality compares identical balls;
    use ``is_congruent`` to compare known digits. Congruence is not equality.
    """
    field: PadicField
    unit: int
    v: int
    precision: int

    def __new__(cls, *args, **kwargs):
        raise TypeError("construct p-adic balls through Qp(p, precision)(value)")

    @classmethod
    def _from_fraction(cls, field, x: Fraction, n: int):
        if type(n) is not int:
            raise TypeError("absolute precision must be an integer")
        obj = object.__new__(cls)
        v = valuation(x, field.p) if x else n
        if v >= n:
            u, v = 0, n
        else:
            scaled = x / p_power(field.p, v)
            modulus = field.p**(n - v)
            u = scaled.numerator * pow(scaled.denominator, -1, modulus) % modulus
        for key, value in (("field", field), ("unit", u), ("v", v), ("precision", n)):
            object.__setattr__(obj, key, value)
        return obj

    @property
    def p(self) -> int:
        """The residue characteristic."""
        return self.field.p

    @property
    def relative_precision(self) -> int:
        """Number of known significant digits (zero for a zero ball)."""
        return self.precision - self.v

    def representative(self) -> Fraction:
        """Return one rational centre, not an inverse embedding into Q or C."""
        return self.unit * p_power(self.p, self.v)

    def valuation(self) -> int:
        """Exact valuation if known; zero balls only give a lower bound."""
        if not self.unit:
            raise PrecisionError(f"only valuation >= {self.precision} is known")
        return self.v

    def norm(self) -> Fraction:
        """Exact p-adic norm for a ball with known nonzero valuation."""
        return p_power(self.p, -self.valuation())

    def is_zero(self) -> bool:
        """Whether zero is consistent with all known digits, not exact equality."""
        return self.unit == 0

    def contains(self, rational) -> bool:
        """Whether an exact rational belongs to this ball."""
        return valuation(self.representative() - exact(rational), self.p) >= self.precision

    def is_congruent(self, other, precision=None) -> bool:
        """Compare modulo p**precision (default: common known precision).

        Raises if either ball lacks the requested precision. An exact rational
        comparison uses all the known digits of this ball.
        """
        if isinstance(other, Padic):
            self._check_prime(other)
            known = min(self.precision, other.precision)
            centre = other.representative()
        else:
            known, centre = self.precision, exact(other)
        n = known if precision is None else precision
        if type(n) is not int:
            raise TypeError("precision must be an integer")
        if n > known:
            raise PrecisionError("requested congruence exceeds known precision")
        return valuation(self.representative() - centre, self.p) >= n

    def truncate(self, precision: int) -> Padic:
        """Forget digits; never increase precision."""
        return self.field(self, absprec=precision)

    def residue(self, exponent=None) -> int:
        """Integer residue modulo p**exponent, for an integral ball."""
        n = self.precision if exponent is None else exponent
        if type(n) is not int or n < 0:
            raise ValueError("exponent must be a nonnegative integer")
        if n > self.precision:
            raise PrecisionError("residue needs additional digits")
        if self.v < 0:
            raise DomainError("a nonintegral p-adic has no integer residue")
        return int(self.representative()) % self.p**n

    def digits(self) -> tuple[int, tuple[int, ...]]:
        """Return (starting exponent, digits in ascending powers of p)."""
        u, out = self.unit, []
        for _ in range(self.relative_precision):
            out.append(u % self.p)
            u //= self.p
        return self.v, tuple(out)

    def _check_prime(self, other):
        if self.p != other.p:
            raise TypeError("arithmetic between different primes is undefined")

    def __add__(self, other):
        if isinstance(other, Padic):
            self._check_prime(other)
            centre, n = other.representative(), min(self.precision, other.precision)
        else:
            try:
                centre, n = exact(other), self.precision
            except TypeError:
                return NotImplemented
        return Padic._from_fraction(self.field, self.representative() + centre, n)

    __radd__ = __add__

    def __neg__(self):
        return Padic._from_fraction(self.field, -self.representative(), self.precision)

    def __sub__(self, other):
        try:
            return self + (-other if isinstance(other, Padic) else -exact(other))
        except TypeError:
            return NotImplemented

    def __rsub__(self, other):
        return -self + other

    def __mul__(self, other):
        if isinstance(other, Padic):
            self._check_prime(other)
            centre = other.representative()
            n = min(self.precision + other.v, other.precision + self.v)
        else:
            try:
                centre = exact(other)
            except TypeError:
                return NotImplemented
            # A finite precision zero is a safe enclosure of exact zero.
            n = self.precision + (valuation(centre, self.p) if centre else 0)
        return Padic._from_fraction(self.field, self.representative() * centre, n)

    __rmul__ = __mul__

    def inverse(self) -> Padic:
        """Invert a ball excluding zero, with precision N - 2*v."""
        if self.is_zero():
            raise PrecisionError("cannot invert a ball containing zero")
        return Padic._from_fraction(self.field, 1 / self.representative(),
                                    self.precision - 2*self.v)

    def __truediv__(self, other):
        if isinstance(other, Padic):
            self._check_prime(other)
            return self * other.inverse()
        try:
            q = exact(other)
        except TypeError:
            return NotImplemented
        if not q:
            raise ZeroDivisionError("division by exact zero")
        return Padic._from_fraction(self.field, self.representative() / q,
                                    self.precision - valuation(q, self.p))

    def __rtruediv__(self, other):
        return self.inverse() * other

    def __pow__(self, exponent: int):
        if type(exponent) is not int:
            raise TypeError("use analytic.power for a principal-unit p-adic exponent")
        if exponent < 0:
            return self.inverse() ** (-exponent)
        if exponent == 0:
            return self.field(1)
        result, base, n = None, self, exponent
        while n:
            if n & 1:
                result = base if result is None else result * base
            n >>= 1
            if n:
                base = base * base
        return result

    def __bool__(self):
        raise TypeError("p-adic truth is ambiguous; use is_zero() or contains()")

    def __eq__(self, other):
        if not isinstance(other, Padic):
            return NotImplemented
        return (self.p, self.unit, self.v, self.precision) == (
            other.p, other.unit, other.v, other.precision)

    def __hash__(self):
        return hash((self.p, self.unit, self.v, self.precision))

    def __repr__(self):
        return f"{self.representative()} + O({self.p}^{self.precision})"
