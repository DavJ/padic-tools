"""Exact polynomials and truncated formal series over Q; optional SymPy export."""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import factorial, isqrt

from .core import DomainError, Padic, PadicField, exact


@dataclass(frozen=True, init=False)
class Polynomial:
    """Univariate Q[x], coefficients in ascending degree order."""
    coefficients: tuple[Fraction, ...]

    def __init__(self, coefficients=()):
        cs = [exact(c) for c in coefficients]
        while cs and not cs[-1]:
            cs.pop()
        object.__setattr__(self, "coefficients", tuple(cs))

    @property
    def degree(self):
        """Degree, using -1 for the zero polynomial."""
        return len(self.coefficients)-1

    def __call__(self, x):
        if not self.coefficients:
            return x.field.O(x.precision) if isinstance(x, Padic) else Fraction(0)
        result = self.coefficients[-1]
        for c in reversed(self.coefficients[:-1]):
            result = result*x+c
        return x.field(result) if isinstance(x, Padic) and not isinstance(result, Padic) else result

    def __add__(self, other):
        other = other if isinstance(other, Polynomial) else Polynomial([other])
        n = max(len(self.coefficients), len(other.coefficients))
        return Polynomial([self[i]+other[i] for i in range(n)])

    __radd__ = __add__

    def __getitem__(self, i):
        return self.coefficients[i] if 0 <= i < len(self.coefficients) else Fraction(0)

    def __neg__(self):
        return Polynomial([-c for c in self.coefficients])

    def __sub__(self, other):
        return self + (-other if isinstance(other, Polynomial) else -exact(other))

    def __rsub__(self, other):
        return -self+other

    def __mul__(self, other):
        other = other if isinstance(other, Polynomial) else Polynomial([other])
        if self.degree < 0 or other.degree < 0:
            return Polynomial()
        cs = [Fraction(0)]*(self.degree+other.degree+1)
        for i, a in enumerate(self.coefficients):
            for j, b in enumerate(other.coefficients):
                cs[i+j] += a*b
        return Polynomial(cs)

    __rmul__ = __mul__

    def __pow__(self, n):
        if type(n) is not int or n < 0:
            raise ValueError("polynomial exponent must be a nonnegative integer")
        result, base = Polynomial([1]), self
        while n:
            if n & 1:
                result = result*base
            n >>= 1
            if n:
                base = base*base
        return result

    def __divmod__(self, divisor):
        divisor = divisor if isinstance(divisor, Polynomial) else Polynomial([divisor])
        if divisor.degree < 0:
            raise ZeroDivisionError("zero polynomial divisor")
        r, q = self, [Fraction(0)]*max(0, self.degree-divisor.degree+1)
        while r.degree >= divisor.degree:
            i, c = r.degree-divisor.degree, r.coefficients[-1]/divisor.coefficients[-1]
            q[i] = c
            r -= Polynomial([0]*i+[c])*divisor
        return Polynomial(q), r

    def derivative(self):
        """Formal derivative in Q[x]."""
        return Polynomial([i*c for i, c in enumerate(self.coefficients) if i])

    def integral(self, constant=0):
        """Formal antiderivative, with an exact integration constant."""
        return Polynomial([constant]+[c/(i+1) for i, c in enumerate(self.coefficients)])

    def gcd(self, other):
        """Monic polynomial gcd via the Euclidean algorithm."""
        a, b = self, other
        while b.degree >= 0:
            a, b = b, divmod(a, b)[1]
        return a*(1/a.coefficients[-1]) if a.degree >= 0 else a

    def to_sympy(self, variable="x"):
        """Convert to an exact SymPy expression; requires the symbolic extra."""
        import sympy as sp
        x = sp.Symbol(variable)
        return sum((sp.Rational(c.numerator, c.denominator)*x**i
                    for i, c in enumerate(self.coefficients)), sp.Integer(0))


@dataclass(frozen=True, init=False)
class PowerSeries:
    """An element of Q[[x]]/(x**order).

    Coefficients beyond ``order`` are unknown, not zero. ``polynomial()``
    explicitly chooses the finite truncation. Formal arithmetic carries
    x-adic order independently of the p-adic precision of numerical values.
    """
    coefficients: tuple[Fraction, ...]
    order: int

    def __init__(self, coefficients=(), order=10):
        if type(order) is not int or order < 1:
            raise ValueError("series order must be a positive integer")
        cs = tuple(exact(c) for c in coefficients)[:order]
        object.__setattr__(self, "coefficients", cs+(Fraction(0),)*(order-len(cs)))
        object.__setattr__(self, "order", order)

    def _coerce(self, other):
        return other if isinstance(other, PowerSeries) else PowerSeries([other], self.order)

    def __add__(self, other):
        other = self._coerce(other)
        n = min(self.order, other.order)
        return PowerSeries([self.coefficients[i]+other.coefficients[i] for i in range(n)], n)

    __radd__ = __add__

    def __neg__(self):
        return PowerSeries([-c for c in self.coefficients], self.order)

    def __sub__(self, other):
        return self + (-self._coerce(other))

    def __rsub__(self, other):
        return -self+other

    def __mul__(self, other):
        other = self._coerce(other)
        n = min(self.order, other.order)
        return PowerSeries([sum((self.coefficients[j]*other.coefficients[i-j]
                                 for j in range(i+1)), Fraction(0)) for i in range(n)], n)

    __rmul__ = __mul__

    def inverse(self):
        """Formal inverse when the constant coefficient is nonzero."""
        a = self.coefficients
        if not a[0]:
            raise ZeroDivisionError("series inverse requires a nonzero constant")
        b = [1/a[0]]
        for n in range(1, self.order):
            b.append(-sum((a[k]*b[n-k] for k in range(1, n+1)), Fraction(0))/a[0])
        return PowerSeries(b, self.order)

    def __truediv__(self, other):
        return self * self._coerce(other).inverse()

    def __pow__(self, n):
        if type(n) is not int:
            raise TypeError("formal series exponent must be an integer")
        if n < 0:
            return self.inverse()**(-n)
        result, base = PowerSeries([1], self.order), self
        while n:
            if n & 1:
                result = result*base
            n >>= 1
            if n:
                base = base*base
        return result

    def derivative(self):
        """Derivative; one coefficient of order is lost."""
        if self.order == 1:
            raise DomainError("derivative has no known coefficients at order 1")
        return PowerSeries([i*c for i, c in enumerate(self.coefficients) if i], self.order-1)

    def integral(self, constant=0):
        """Formal integration increases the known order by one."""
        return PowerSeries([constant]+[c/(i+1) for i, c in enumerate(self.coefficients)], self.order+1)

    def exp(self):
        """Formal exponential, requiring zero constant over Q."""
        a = self.coefficients
        if a[0]:
            raise DomainError("formal exp over Q requires constant coefficient zero")
        b = [Fraction(1)]
        for n in range(1, self.order):
            b.append(sum((k*a[k]*b[n-k] for k in range(1, n+1)), Fraction(0))/n)
        return PowerSeries(b, self.order)

    def log(self):
        """Formal logarithm, requiring constant coefficient one."""
        if self.coefficients[0] != 1:
            raise DomainError("formal log over Q requires constant coefficient one")
        if self.order == 1:
            return PowerSeries([], 1)
        return (self.derivative()*self.inverse()).integral()

    def compose(self, inner):
        """Compose formal series; inner must have zero constant coefficient."""
        if not isinstance(inner, PowerSeries) or inner.coefficients[0]:
            raise DomainError("formal composition requires an inner series with constant zero")
        n = min(self.order, inner.order)
        result = PowerSeries([], n)
        for c in reversed(self.coefficients[:n]):
            result = result*inner+c
        return result

    def polynomial(self) -> Polynomial:
        """The chosen finite polynomial truncation; no analytic tail is asserted."""
        return Polynomial(self.coefficients)

    def to_sympy(self, variable="x"):
        """Exact SymPy series including its formal O(x**order) term."""
        import sympy as sp
        x = sp.Symbol(variable)
        return self.polynomial().to_sympy(variable)+sp.Order(x**self.order, (x, 0))


def exp_series(order=10):
    """Exact formal Taylor expansion of exp(x) through order-1."""
    return PowerSeries([Fraction(1, factorial(k)) for k in range(order)], order)


def log1p_series(order=10):
    """Exact formal Taylor expansion of log(1+x)."""
    return PowerSeries([0]+[Fraction((-1)**(k+1), k) for k in range(1, order)], order)


def theta3_series(order=20):
    """Formal theta constant 1+2*q+2*q**4+..., with q-adic order."""
    cs = [Fraction(0)]*order
    if order < 1:
        raise ValueError("order must be positive")
    cs[0] = Fraction(1)
    for k in range(1, isqrt(order-1)+1):
        cs[k*k] = Fraction(2)
    return PowerSeries(cs, order)


def hensel_lift(f: Polynomial, root: int, field: PadicField, *, precision=None) -> Padic:
    """Lift a simple root modulo p to absolute precision N.

    Requires p-integral coefficients and f'(root) nonzero modulo p.
    Singular roots are rejected instead of assuming the simple-root theorem.
    """
    from .core import valuation
    n = field.precision if precision is None else precision
    if type(n) is not int or n < 1 or type(root) is not int:
        raise ValueError("root must be an integer and precision positive")
    p, r, df = field.p, root % field.p, f.derivative()
    if any(valuation(c, p) < 0 for c in f.coefficients):
        raise DomainError("Hensel lift requires p-integral coefficients")
    def mod(q, m):
        return q.numerator*pow(q.denominator, -1, m) % m
    if mod(f(Fraction(r)), p) != 0 or mod(df(Fraction(r)), p) == 0:
        raise DomainError("requires f(root)=0 and f'(root)!=0 modulo p")
    modulus = p
    for _ in range(1, n):
        next_modulus = modulus*p
        digit = -(mod(f(Fraction(r)), next_modulus)//modulus)*pow(mod(df(Fraction(r)), p), -1, p) % p
        r += digit*modulus
        modulus = next_modulus
    return field(r, absprec=n)
