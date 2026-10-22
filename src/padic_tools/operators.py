"""Whole-Qp Vladimirov derivatives and conservative ball diffusion.

Positive integer orders retain exact rational spectral/kernel constants.
The whole-space image includes its analytic exterior tail. The regional
ball generator integrates jumps only within the chosen compact ball.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from ._cyclotomic import CyclotomicValue, _size
from .core import DomainError, Padic, exact, p_power, valuation
from .harmonic import PadicTestFunction, finite_fourier, haar_volume


def _order(order):
    if type(order) is not int or order < 1:
        raise ValueError("Vladimirov order must be a positive integer")
    return order


@dataclass(frozen=True, init=False)
class BallVladimirov:
    """Conservative regional generator on B=p**a*Z_p modulo p**b*Z_p.

    ``Lf=C*integral_B (f(x)-f(y))/abs_p(x-y)**(1+order) dy``,
    C=(p**order-1)/(1-p**(-1-order)). Its constant mode has eigenvalue zero.
    Other frequencies have eigenvalue ``abs_p(xi)**order-escape_rate``. This is
    distinct from the whole-Qp derivative of a zero-extended function.
    """
    template: PadicTestFunction
    order: int
    kernel_constant: Fraction
    escape_rate: Fraction
    eigenvalues: tuple[Fraction, ...]

    def __init__(self, field, order=1, *, support_exponent=0,
                 constancy_exponent=1, max_size=256):
        order = _order(order)
        template = PadicTestFunction.ball(field)
        p,a,b = template.field.p,support_exponent,constancy_exponent
        if type(a) is not int or type(b) is not int:
            raise TypeError("lattice exponents must be integers")
        if b < a:
            raise ValueError("constancy exponent must be >= support exponent")
        size = _size(p,b-a,max_size)
        template = PadicTestFunction(template.field,[0]*size,support_exponent=a,
                                    constancy_exponent=b,max_size=max_size)
        constant = (p**order-1)/(1-p_power(p,-1-order))
        escape = (1-Fraction(1,p))*p_power(p,order*a)/(1-p_power(p,-1-order))
        eigenvalues = (Fraction(0),)+tuple(p_power(p,order*(b-valuation(k,p)))-escape
                                          for k in range(1,size))
        for name,value in (("template",template),("order",order),
                           ("kernel_constant",constant),("escape_rate",escape),
                           ("eigenvalues",eigenvalues)):
            object.__setattr__(self,name,value)

    def apply(self, function):
        """Apply the exact regional generator via its finite Fourier spectrum."""
        self.template._check(function)
        size = len(function.values)
        spectrum = finite_fourier(function.values,function.field,max_size=size)
        values = finite_fourier([value*rate for value,rate in zip(spectrum,self.eigenvalues)],
                                function.field,inverse=True,max_size=size)
        return function._replace_values(values)

    def matrix(self):
        """Independent exact kernel matrix: rows sum to zero, off-diagonals <0.

        -matrix is the finite-state continuous-time Markov generator.
        Each off-diagonal integral includes the real Haar coset weight.
        """
        p = self.template.field.p
        a,b = self.template.support_exponent,self.template.constancy_exponent
        size = len(self.template.values)
        weight = self.kernel_constant*haar_volume(p,b)
        rows = []
        for j in range(size):
            row = [-weight*p_power(p,(1+self.order)*(a+valuation(j-k,p)))
                   if j != k else Fraction(0) for k in range(size)]
            row[j] = -sum(row)
            rows.append(tuple(row))
        return tuple(rows)

    def evolve(self, function, time, *, dps=50):
        """Numerical exp(-time*L) using isolated mpmath (complex extra).

        Time is exact nonnegative real rational data. Working precision is
        dps; outputs are archimedean approximations, not certified balls.
        The exact constant mode is retained by the mathematical semigroup.
        """
        from .bridges import _complex_context, _mp_rational
        self.template._check(function)
        time = exact(time)
        if time < 0:
            raise DomainError("diffusion time must be nonnegative")
        output = _complex_context(dps)
        ctx = _complex_context(dps+10)
        size = len(function.values)
        spectrum = finite_fourier(function.values,function.field,max_size=size)
        t = _mp_rational(ctx,time)
        decayed = [value.complex(dps+10)*ctx.exp(-t*_mp_rational(ctx,rate))
                   for value,rate in zip(spectrum,self.eigenvalues)]
        values = tuple(output.mpc(ctx.fsum(value*ctx.expjpi(ctx.mpf(2*j*k)/size)
                                          for k,value in enumerate(decayed))) for j in range(size))
        return DiffusionState(function,time,dps,values)


@dataclass(frozen=True)
class DiffusionState:
    """Approximate complex values of a conservative finite-ball heat solution.

    template records the exact initial data and spatial lattice. dps is
    archimedean working precision; it is unrelated to field.precision.
    """
    template: PadicTestFunction
    time: Fraction
    dps: int
    values: tuple

    def __call__(self, point):
        """Sample on the declared lattice, requiring sufficient point digits."""
        index = self.template._index(point)
        return 0 if index is None else self.values[index]

    def integral(self):
        """Approximate real-Haar integral of the complex output."""
        from .bridges import _complex_context, _mp_rational
        ctx = _complex_context(self.dps)
        return ctx.fsum(self.values)*_mp_rational(ctx,haar_volume(
            self.template.field,self.template.constancy_exponent))

    def norm_squared(self):
        """Approximate Haar L2 norm squared."""
        from .bridges import _complex_context, _mp_rational
        ctx = _complex_context(self.dps)
        return ctx.fsum(ctx.fabs(value)**2 for value in self.values)*_mp_rational(
            ctx,haar_volume(self.template.field,self.template.constancy_exponent))


@dataclass(frozen=True)
class VladimirovImage:
    """Exact whole-Qp derivative: finite interior plus radial exterior tail.

    Outside ``B=p**a*Z_p``, the value is ``tail_coefficient*abs_p(x)**(-1-order)``.
    This tail is part of the result, including when the initial function
    has nonzero mass. Zero-mean inputs have a compactly supported image.
    """
    interior: PadicTestFunction
    order: int
    tail_coefficient: CyclotomicValue

    def __call__(self, point):
        index = self.interior._index(point)
        if index is not None:
            return self.interior.values[index]
        if isinstance(point,Padic):
            v = point.v
        else:
            v = valuation(point,self.interior.field.p)
        return self.tail_coefficient*p_power(self.interior.field.p,(1+self.order)*v)

    def as_test_function(self):
        """Return the compact image only if its exact exterior tail is zero."""
        if self.tail_coefficient:
            raise DomainError("whole-Qp Vladimirov image has a noncompact exterior tail")
        return self.interior

    def integral(self):
        """Exact integral including all exterior shells."""
        p,a = self.interior.field.p,self.interior.support_exponent
        shell_sum = (1-Fraction(1,p))*p_power(p,self.order*(a-1))/(1-p_power(p,-self.order))
        return self.interior.integral()+self.tail_coefficient*shell_sum

    def norm_squared(self):
        """Exact L2 norm including the full square-integrable exterior tail."""
        p,a = self.interior.field.p,self.interior.support_exponent
        exponent = 1+2*self.order
        shell_sum = (1-Fraction(1,p))*p_power(p,exponent*(a-1))/(1-p_power(p,-exponent))
        return self.interior.norm_squared()+self.tail_coefficient*self.tail_coefficient.conjugate()*shell_sum


def vladimirov(function, order=1, *, max_size=256):
    """Exact Fourier multiplier abs_p(xi)**order on Qp, including its full tail.

    Positive integer orders and exact PadicTestFunction source values are
    supported. Inside the source container, the result is regional Lf plus
    ``escape_rate*f``. Outside it equals ``-C*integral(f)*abs_p(x)**(-1-order)``.
    """
    if not isinstance(function,PadicTestFunction):
        raise TypeError("expected a PadicTestFunction")
    operator = BallVladimirov(function.field,order,
        support_exponent=function.support_exponent,
        constancy_exponent=function.constancy_exponent,max_size=max_size)
    interior = operator.apply(function)+function*operator.escape_rate
    return VladimirovImage(interior,operator.order,-function.integral()*operator.kernel_constant)
