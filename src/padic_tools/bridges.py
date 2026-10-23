"""Explicit diagonal embeddings and scalar base changes.

There is no continuous field bridge Q_p -> C fixing Q: p**n tends to
zero p-adically and grows without bound in C. We instead map a common
exact rational algebra separately into the two coefficient domains.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from .core import DomainError, PadicField, exact, valuation
from .analytic import sqrt
from .algebra import Gaussian, QuadraticAlgebra, QuadraticElement, Quaternion, QuaternionAlgebra
from .extensions import UnramifiedQuadraticField
from .ramified import EisensteinQuadraticField


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


def _local_target(target, absprec):
    field = target.field if isinstance(target, GaussianLocalPlace) else target
    if not isinstance(field, (PadicField, UnramifiedQuadraticField, EisensteinQuadraticField)):
        raise TypeError("target must be a certified local field or GaussianLocalPlace")
    n = field.precision if absprec is None else absprec
    if type(n) is not int:
        raise TypeError("absprec must be an integer in target uniformizer units")
    return field, n


@dataclass(frozen=True, init=False)
class GaussianLocalPlace:
    """A chosen embedding of exact Q[i] into its completion at p.

    ``root_sign`` is required and is either +1 or -1. At split primes the
    positive root lifts the smaller residue of sqrt(-1) modulo p. At inert
    primes it is the recorded quadratic generator. At p=2 it is pi-1,
    where pi**2-2*pi+2=0. The other sign gives the conjugate embedding.

    ``absprec`` counts absolute powers of the target uniformizer: p for
    odd primes, pi for p=2. Exact sources can be reevaluated at any precision.
    """
    base: PadicField
    root_sign: int
    field: object

    def __init__(self, base: PadicField, *, root_sign):
        if not isinstance(base, PadicField):
            raise TypeError("base must be Qp")
        if type(root_sign) is not int or root_sign not in (-1, 1):
            raise ValueError("root_sign must be the integer +1 or -1")
        if base.p == 2:
            field = EisensteinQuadraticField(base, [2, -2, 1], name="pi")
        elif base.p % 4 == 3:
            field = UnramifiedQuadraticField(base, [1, 0, 1], name="i")
        else:
            field = base
        object.__setattr__(self, "base", base)
        object.__setattr__(self, "root_sign", root_sign)
        object.__setattr__(self, "field", field)

    @property
    def kind(self):
        """One of 'split', 'unramified', 'ramified'."""
        return "ramified" if self.base.p == 2 else ("split" if self.base.p % 4 == 1 else "unramified")

    @property
    def ramification_index(self):
        """The target's ramification index over Qp."""
        return 2 if self.base.p == 2 else 1

    @property
    def residue_degree(self):
        """The target's residue degree over Fp."""
        return 2 if self.base.p % 4 == 3 else 1

    def __call__(self, value=0, imag=0, *, absprec=None):
        """Evaluate an exact Gaussian source, rounding once in the target."""
        z = GaussianBridge(value, imag).value
        field, n = _local_target(self, absprec)
        a, b = z.a, self.root_sign*z.b
        if self.base.p == 2:
            return field(a-b, b, absprec=n)
        if self.kind == "unramified":
            return field(a, b, absprec=n)
        if not b:
            return field(a, absprec=n)
        # b can lose digits; choose root precision before multiplying by b.
        root_precision = max(1, n-valuation(b, self.base.p))
        root = sqrt(self.base(-1, absprec=root_precision))
        return field(a+b*root, absprec=n)

    def generator_image(self, *, absprec=None):
        """The chosen image of the central Gaussian i."""
        return self(0, 1, absprec=absprec)

    def conjugate(self):
        """The opposite embedding; at split primes this swaps the two places."""
        return GaussianLocalPlace(self.base, root_sign=-self.root_sign)


@dataclass(frozen=True, init=False)
class RationalBridge:
    """One exact q in Q, with independent complex and p-adic images."""
    value: Fraction

    def __init__(self, value):
        object.__setattr__(self, "value", exact(value))

    def padic(self, field: PadicField):
        """The image in the specified finite precision Q_p."""
        return field(self.value)

    def local(self, target, *, absprec=None):
        """Image in a certified local field, at absolute uniformizer precision."""
        field, n = _local_target(target, absprec)
        return field(self.value, absprec=n)

    def complex(self, dps=50):
        """Complex image using an isolated mpmath context (complex extra)."""
        ctx = _complex_context(dps)
        return ctx.mpc(_mp_rational(ctx, self.value))


@dataclass(frozen=True, init=False)
class GaussianBridge:
    """An exact element of Q[i], with explicit choices of local image.

    ``padic_algebra`` always returns the quotient algebra Q_p[t]/(t^2+1).
    ``split_padic`` chooses a root of -1 in Q_p, possible only for p=1 mod 4.
    ``local`` uses a GaussianLocalPlace and returns an actual local-field
    ball, including ramified Q2(i), with uniformizer precision.
    """
    value: QuadraticElement

    def __init__(self, value=0, imag=0):
        if isinstance(value, QuadraticElement):
            if value.parent != Gaussian:
                raise TypeError("source must be exact Q[i], not p-adic coefficients")
            if exact(imag):
                raise TypeError("imag must be zero when passing a Gaussian element")
            z = Gaussian(exact(value.a), exact(value.b))
        else:
            z = Gaussian(value, imag)
        object.__setattr__(self, "value", z)

    def local(self, place: GaussianLocalPlace, *, absprec=None):
        """Evaluate in the specified place, with absolute target precision."""
        if not isinstance(place, GaussianLocalPlace):
            raise TypeError("Gaussian images require an explicit GaussianLocalPlace")
        return place(self.value, absprec=absprec)

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
        if (not isinstance(self.value, Quaternion)
                or not isinstance(self.value.parent, QuaternionAlgebra)
                or self.value.parent.base not in (exact, Gaussian)):
            raise TypeError("source must be a quaternion over Q or exact Q[i]")
        if len(self.value.coefficients) != 4:
            raise TypeError("a quaternion source requires four exact coefficients")
        converter = exact if self.value.parent.base == exact else lambda c: GaussianBridge(c).value
        coefficients = tuple(converter(c) for c in self.value.coefficients)
        object.__setattr__(self, "value", Quaternion(self.value.parent, coefficients))

    def local(self, target, *, absprec=None):
        """Base change to a certified field, choosing a place for Q[i].

        All four coefficients start at the same absolute target precision.
        Quaternion conjugation and reduced norm stay over that scalar field.
        """
        field, n = _local_target(target, absprec)
        if self.value.parent.base == Gaussian:
            if not isinstance(target, GaussianLocalPlace):
                raise TypeError("biquaternion coefficients require a GaussianLocalPlace")
            mapper = lambda c: target(c, absprec=n)
        else:
            mapper = lambda c: field(c, absprec=n)
        algebra = QuaternionAlgebra(field, self.value.parent.a, self.value.parent.b)
        return self.value.map_coefficients(algebra, mapper)

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
