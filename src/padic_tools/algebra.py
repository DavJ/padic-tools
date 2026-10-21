"""Quadratic, quaternion and diagonal Clifford algebras over a common base.

These are algebraic constructions. Replacing the coefficient field does not
identify its topology with that of any other completion.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from types import MappingProxyType

from .core import DomainError, Padic, exact


def _exact_zero(x):
    # O(p^N) is an uncertainty, so it must not be pruned as exact zero.
    if isinstance(x, Padic):
        return False
    if isinstance(x, QuadraticElement):
        return _exact_zero(x.a) and _exact_zero(x.b)
    return x == 0


def _congruent(a, b, precision=None):
    return a.is_congruent(b, precision) if hasattr(a, "is_congruent") else a == b


@dataclass(frozen=True)
class QuadraticAlgebra:
    """The commutative quotient base[t]/(t**2-d).

    This need not be a field: at a split prime Q_p[i] is Q_p x Q_p.
    Inversion checks the norm. ``base`` is a scalar conversion callable.
    """
    base: object = exact
    d: Fraction = Fraction(-1)
    name: str = "i"

    def __post_init__(self):
        object.__setattr__(self, "d", exact(self.d))
        if not self.d:
            raise DomainError("quadratic parameter d must be nonzero")

    def __call__(self, a=0, b=0):
        if isinstance(a, QuadraticElement):
            if a.parent == self and _exact_zero(b):
                return a
            raise TypeError("use map_coefficients for different quadratic algebras")
        return QuadraticElement(self, self.base(a), self.base(b))

    def generator(self):
        """The chosen formal square root of d."""
        return self(0, 1)


@dataclass(frozen=True)
class QuadraticElement:
    """a+b*t, with t**2=d and a recorded parent algebra."""
    parent: QuadraticAlgebra
    a: object
    b: object

    def _coerce(self, other):
        if isinstance(other, QuadraticElement):
            if self.parent != other.parent:
                raise TypeError("quadratic parents must agree")
            return other
        return self.parent(other)

    def __add__(self, other):
        try:
            other = self._coerce(other)
        except TypeError:
            return NotImplemented
        return QuadraticElement(self.parent, self.a+other.a, self.b+other.b)

    __radd__ = __add__

    def __neg__(self):
        return QuadraticElement(self.parent, -self.a, -self.b)

    def __sub__(self, other):
        return self+-self._coerce(other)

    def __rsub__(self, other):
        return -self+other

    def __mul__(self, other):
        try:
            other = self._coerce(other)
        except TypeError:
            return NotImplemented
        return QuadraticElement(self.parent,
            self.a*other.a+self.parent.d*self.b*other.b,
            self.a*other.b+self.b*other.a)

    __rmul__ = __mul__

    def conjugate(self):
        """Quadratic conjugation a+b*t -> a-b*t."""
        return QuadraticElement(self.parent, self.a, -self.b)

    def norm(self):
        """Norm a**2-d*b**2 in the base algebra."""
        return self.a*self.a-self.parent.d*self.b*self.b

    def inverse(self):
        """Invert if the norm is a known invertible base element."""
        n = self.norm()
        if _exact_zero(n):
            raise ZeroDivisionError("quadratic element has zero norm")
        return QuadraticElement(self.parent, self.a/n, -self.b/n)

    def __truediv__(self, other):
        return self*self._coerce(other).inverse()

    def __rtruediv__(self, other):
        return self.inverse()*other

    def __pow__(self, n):
        if type(n) is not int:
            raise TypeError("quadratic exponent must be an integer")
        if n < 0:
            return self.inverse()**(-n)
        result, base = self.parent(1), self
        while n:
            if n & 1:
                result = result*base
            n >>= 1
            if n:
                base = base*base
        return result

    def map_coefficients(self, target: QuadraticAlgebra, mapper=None):
        """Change base via a specified coefficient map, preserving t**2=d."""
        if self.parent.d != target.d:
            raise DomainError("coefficient maps must preserve the defining polynomial")
        f = target.base if mapper is None else mapper
        return QuadraticElement(target, f(self.a), f(self.b))

    def is_congruent(self, other, precision=None):
        other = self._coerce(other)
        return _congruent(self.a, other.a, precision) and _congruent(self.b, other.b, precision)

    def __repr__(self):
        return f"({self.a}) + ({self.b})*{self.parent.name}"


Gaussian = QuadraticAlgebra()


@dataclass(frozen=True)
class QuaternionAlgebra:
    """Quaternion algebra (a,b) over a commutative characteristic-zero base.

    i**2=a, j**2=b, i*j=k=-j*i. Hamilton parameters are a=b=-1.
    A p-adic Hamilton algebra can be split; a nonzero element need not invert.
    """
    base: object = exact
    a: Fraction = Fraction(-1)
    b: Fraction = Fraction(-1)

    def __post_init__(self):
        object.__setattr__(self, "a", exact(self.a))
        object.__setattr__(self, "b", exact(self.b))
        if not self.a or not self.b:
            raise DomainError("quaternion parameters must be nonzero")

    def __call__(self, w=0, x=0, y=0, z=0):
        if isinstance(w, Quaternion):
            if w.parent == self and all(_exact_zero(c) for c in (x, y, z)):
                return w
            raise TypeError("use map_coefficients for different quaternion algebras")
        return Quaternion(self, tuple(self.base(c) for c in (w, x, y, z)))

    def generators(self):
        """The three chosen quaternion units (i,j,k)."""
        return self(0, 1), self(0, 0, 1), self(0, 0, 0, 1)


@dataclass(frozen=True)
class Quaternion:
    """w+x*i+y*j+z*k over a recorded quaternion algebra."""
    parent: QuaternionAlgebra
    coefficients: tuple

    def _coerce(self, other):
        if isinstance(other, Quaternion):
            if self.parent != other.parent:
                raise TypeError("quaternion parents must agree")
            return other
        return self.parent(other)

    def __add__(self, other):
        other = self._coerce(other)
        return Quaternion(self.parent, tuple(a+b for a, b in zip(self.coefficients, other.coefficients)))

    __radd__ = __add__

    def __neg__(self):
        return Quaternion(self.parent, tuple(-a for a in self.coefficients))

    def __sub__(self, other):
        return self+-self._coerce(other)

    def __rsub__(self, other):
        return -self+other

    def __mul__(self, other):
        other = self._coerce(other)
        w, x, y, z = self.coefficients
        W, X, Y, Z = other.coefficients
        a, b = self.parent.a, self.parent.b
        return Quaternion(self.parent, (
            w*W+a*x*X+b*y*Y-a*b*z*Z,
            w*X+x*W+b*(z*Y-y*Z),
            w*Y+y*W+a*(x*Z-z*X),
            w*Z+z*W+x*Y-y*X))

    __rmul__ = __mul__

    def conjugate(self):
        """Quaternion conjugation: negate the three nonscalar components."""
        w, x, y, z = self.coefficients
        return Quaternion(self.parent, (w, -x, -y, -z))

    def norm(self):
        """Reduced norm w**2-a*x**2-b*y**2+a*b*z**2."""
        w, x, y, z = self.coefficients
        a, b = self.parent.a, self.parent.b
        return w*w-a*x*x-b*y*y+a*b*z*z

    def inverse(self):
        """Conjugate divided by reduced norm; zero divisors are rejected."""
        n = self.norm()
        if _exact_zero(n):
            raise ZeroDivisionError("quaternion has zero reduced norm")
        return Quaternion(self.parent, tuple(c/n for c in self.conjugate().coefficients))

    def __truediv__(self, other):
        return self*self._coerce(other).inverse()

    def __rtruediv__(self, other):
        return self.inverse()*other

    def __pow__(self, n):
        if type(n) is not int:
            raise TypeError("quaternion exponent must be an integer")
        if n < 0:
            return self.inverse()**(-n)
        result, base = self.parent(1), self
        while n:
            if n & 1:
                result = result*base
            n >>= 1
            if n:
                base = base*base
        return result

    def map_coefficients(self, target: QuaternionAlgebra, mapper=None):
        """Base change preserving the quaternion defining relations."""
        if (self.parent.a, self.parent.b) != (target.a, target.b):
            raise DomainError("quaternion parameters must be preserved")
        f = target.base if mapper is None else mapper
        return Quaternion(target, tuple(f(c) for c in self.coefficients))

    def is_congruent(self, other, precision=None):
        other = self._coerce(other)
        return all(_congruent(a, b, precision) for a, b in zip(self.coefficients, other.coefficients))

    def to_clifford(self):
        """Hamilton algebra -> Cl(3,0)^even, over the same coefficient base.

        i=e2*e3, j=e3*e1, k=-e1*e2. Complex coefficients give the same
        identification for biquaternions. This is an algebra isomorphism.
        """
        if (self.parent.a, self.parent.b) != (-1, -1):
            raise DomainError("this explicit Clifford map requires Hamilton parameters")
        C = CliffordAlgebra((1, 1, 1), self.parent.base)
        w, x, y, z = self.coefficients
        return C({0: w, 6: x, 5: -y, 3: -z})


Hamilton = QuaternionAlgebra()
Biquaternions = QuaternionAlgebra(Gaussian)


@dataclass(frozen=True, init=False)
class CliffordAlgebra:
    """Diagonal Clifford algebra, e_i**2=metric[i], e_i*e_j=-e_j*e_i.

    Blade keys are integer bitmasks, e.g. mask 0=scalar, 3=e1*e2.
    The coefficient base can be Q, Q_p, a quadratic algebra or another
    commutative characteristic-zero scalar domain.
    """
    metric: tuple
    base: object

    def __init__(self, metric=(1, 1, 1), base=exact):
        object.__setattr__(self, "metric", tuple(exact(q) for q in metric))
        object.__setattr__(self, "base", base)
        if len(self.metric) > 12:
            raise DomainError("native sparse Clifford backend supports dimension <= 12")

    def __call__(self, coefficients=0):
        if isinstance(coefficients, Multivector):
            if coefficients.parent != self:
                raise TypeError("Clifford parents must agree")
            return coefficients
        cs = coefficients if hasattr(coefficients, "items") else {0: coefficients}
        out = {}
        for mask, c in cs.items():
            if type(mask) is not int or not 0 <= mask < 2**len(self.metric):
                raise ValueError("blade bitmask exceeds Clifford dimension")
            cc = self.base(c)
            if not _exact_zero(cc):
                out[mask] = cc
        return Multivector(self, MappingProxyType(out))

    def e(self, index):
        """Basis vector, using zero-based indices."""
        if type(index) is not int or not 0 <= index < len(self.metric):
            raise IndexError("basis vector index out of range")
        return self({1 << index: 1})

    def pseudoscalar(self):
        """The oriented product e1*...*en."""
        return self({2**len(self.metric)-1: 1})

    def from_quaternion(self, value: Quaternion):
        """Hamilton quaternion in this Cl(3,0) even subalgebra."""
        image = value.to_clifford()
        if image.parent != self:
            raise DomainError("requires Cl(3,0) over the same scalar base")
        return image


@dataclass(frozen=True)
class Multivector:
    """Sparse multivector with geometric and exterior products."""
    parent: CliffordAlgebra
    coefficients: object

    def _coerce(self, other):
        return self.parent(other)

    def __add__(self, other):
        other = self._coerce(other)
        cs = dict(self.coefficients)
        for mask, c in other.coefficients.items():
            cs[mask] = cs[mask]+c if mask in cs else c
        return self.parent(cs)

    __radd__ = __add__

    def __neg__(self):
        return self.parent({mask: -c for mask, c in self.coefficients.items()})

    def __sub__(self, other):
        return self+-self._coerce(other)

    def __rsub__(self, other):
        return -self+other

    def _product(self, other, wedge=False):
        other = self._coerce(other)
        cs, metric = {}, self.parent.metric
        for a, ca in self.coefficients.items():
            for b, cb in other.coefficients.items():
                if wedge and a & b:
                    continue
                sign = (-1)**sum((b & ((1 << i)-1)).bit_count()
                                  for i in range(len(metric)) if a & (1 << i))
                factor = Fraction(sign)
                for i, q in enumerate(metric):
                    if (a & b) & (1 << i):
                        factor *= q
                mask, term = a ^ b, ca*cb*factor
                cs[mask] = cs[mask]+term if mask in cs else term
        return self.parent(cs)

    def __mul__(self, other):
        return self._product(other)

    def __rmul__(self, other):
        return self._coerce(other)._product(self)

    def wedge(self, other):
        """Exterior product, also available as a ^ b."""
        return self._product(other, wedge=True)

    __xor__ = wedge

    def grade(self, k):
        """Projection to grade k."""
        return self.parent({m: c for m, c in self.coefficients.items() if m.bit_count() == k})

    def reverse(self):
        """Reverse basis-vector order in every blade."""
        return self.parent({m: c*(-1)**(m.bit_count()*(m.bit_count()-1)//2)
                            for m, c in self.coefficients.items()})

    def involute(self):
        """Grade involution: negate odd grades."""
        return self.parent({m: c*(-1)**m.bit_count() for m, c in self.coefficients.items()})

    def conjugate(self):
        """Clifford conjugation, the composition of reverse and involution."""
        return self.reverse().involute()

    def is_congruent(self, other, precision=None):
        other = self._coerce(other)
        zero = self.parent.base(0)
        return all(_congruent(self.coefficients.get(m, zero), other.coefficients.get(m, zero), precision)
                   for m in set(self.coefficients) | set(other.coefficients))

    def to_quaternion(self):
        """Inverse of the explicit Hamilton -> Cl(3,0)^even map."""
        if self.parent.metric != (1, 1, 1):
            raise DomainError("requires Cl(3,0)")
        if any(m not in (0, 3, 5, 6) for m in self.coefficients):
            raise DomainError("only the even subalgebra maps to Hamilton quaternions")
        H, c = QuaternionAlgebra(self.parent.base), self.coefficients
        z = self.parent.base(0)
        return H(c.get(0, z), c.get(6, z), -c.get(5, z), -c.get(3, z))
