"""Explicit diagonal embeddings and scalar base changes.

There is no continuous field bridge Q_p -> C fixing Q: p**n tends to
zero p-adically and grows without bound in C. We instead map a common
exact rational algebra separately into the two coefficient domains.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from .core import DomainError, PadicField, exact
from .analytic import sqrt
from .algebra import Gaussian, QuadraticAlgebra, QuadraticElement, Quaternion, QuaternionAlgebra


def _complex_context(dps):
    if type(dps) is not int or dps < 2:
        raise ValueError("dps must be an integer >=2")
    import mpmath
    context = mpmath.mp.clone()
    context.dps = dps
    return context


def _mp_rational(ctx, value):
    q = exact(value)
    return ctx.mpf(q.numerator)/q.denominator


@dataclass(frozen=True, init=False)
class RationalBridge:
    """One exact q in Q, with independent complex and p-adic images."""
    value: Fraction

    def __init__(self, value):
        object.__setattr__(self, "value", exact(value))

    def padic(self, field: PadicField):
        """The image in the specified finite precision Q_p."""
        return field(self.value)

    def complex(self, dps=50):
        """Complex image using an isolated mpmath context (complex extra)."""
        ctx = _complex_context(dps)
        return ctx.mpc(_mp_rational(ctx, self.value))


@dataclass(frozen=True, init=False)
class GaussianBridge:
    """An exact element of Q[i], with two distinct kinds of local image.

    ``padic_algebra`` always returns the quotient algebra Q_p[t]/(t^2+1).
    ``split_padic`` chooses a root of -1 in Q_p, possible only for p=1 mod 4.
    At p=2 the quadratic field is ramified; no uniformizer arithmetic is
    implemented by this coefficient-wise quotient representation.
    """
    value: QuadraticElement

    def __init__(self, value=0, imag=0):
        z = value if isinstance(value, QuadraticElement) else Gaussian(value, imag)
        if z.parent != Gaussian:
            raise TypeError("source must be exact Q[i], not p-adic coefficients")
        object.__setattr__(self, "value", z)

    def complex(self, dps=50):
        """Chosen archimedean embedding sending i to the usual complex i."""
        ctx = _complex_context(dps)
        return ctx.mpc(_mp_rational(ctx, self.value.a), _mp_rational(ctx, self.value.b))

    def padic_algebra(self, field: PadicField):
        """Base change to the quadratic quotient; no false field claim."""
        return self.value.map_coefficients(QuadraticAlgebra(field, -1, "i"))

    def split_padic(self, field: PadicField, *, root_sign=1):
        """Choose i -> +/-sqrt(-1) in Q_p when p=1 mod 4.

        root_sign records the place choice. A root is lifted with guard digits
        to account for a nonintegral imaginary coefficient.
        """
        if field.p % 4 != 1:
            raise DomainError("-1 has no square root in this Q_p; use padic_algebra")
        if root_sign not in (-1, 1) or isinstance(root_sign, bool):
            raise ValueError("root_sign must be +1 or -1")
        from .core import valuation
        b = self.value.b
        guard = max(0, -valuation(b, field.p)) if b else 0
        root = sqrt(field(-1, absprec=field.precision+guard)) * root_sign
        return root*b+self.value.a


@dataclass(frozen=True)
class QuaternionBridge:
    """Base change of exact rational quaternions or Q[i]-biquaternions."""
    value: Quaternion

    def __post_init__(self):
        if not isinstance(self.value, Quaternion) or self.value.parent.base not in (exact, Gaussian):
            raise TypeError("source must be a quaternion over Q or exact Q[i]")

    def padic(self, field: PadicField):
        """Coefficient-wise image over Q_p or Q_p[t]/(t^2+1)."""
        value = self.value
        if value.parent.base == Gaussian:
            target_base = QuadraticAlgebra(field, -1, "i")
            mapper = lambda c: c.map_coefficients(target_base)
        else:
            target_base, mapper = field, field
        target = QuaternionAlgebra(target_base, value.parent.a, value.parent.b)
        return value.map_coefficients(target, mapper)

    def complex(self, dps=50):
        """Coefficient-wise complex quaternion, keeping four components."""
        ctx = _complex_context(dps)
        if self.value.parent.base == Gaussian:
            mapper = lambda c: ctx.mpc(_mp_rational(ctx, c.a), _mp_rational(ctx, c.b))
        else:
            mapper = lambda c: ctx.mpc(_mp_rational(ctx, c))
        target = QuaternionAlgebra(ctx.mpc, self.value.parent.a, self.value.parent.b)
        return self.value.map_coefficients(target, mapper)

    def clifford(self):
        """Exact Hamilton/biquaternion image in Cl(3,0)^even."""
        return self.value.to_clifford()


def riemann_zeta(s, *, dps=50):
    """Classical complex Riemann zeta via mpmath; distinct from zeta_p.

    For exact rational or Gaussian input the common algebraic source is
    preserved up to the requested archimedean rounding. No p-adic ball is
    implicitly converted to a complex value.
    """
    ctx = _complex_context(dps)
    if isinstance(s, QuadraticElement) and s.parent == Gaussian:
        z = ctx.mpc(_mp_rational(ctx, s.a), _mp_rational(ctx, s.b))
    elif isinstance(s, (int, Fraction, str)):
        z = _mp_rational(ctx, s)
    elif isinstance(s, (float, complex)):
        raise TypeError("use rational/Gaussian exact input or a decimal string")
    else:
        raise TypeError("no implicit map from a p-adic or foreign approximation to C")
    return ctx.zeta(z)
