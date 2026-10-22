"""Exact amplitude-one Kozyrev wavelets and finite orthogonal expansions."""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from itertools import islice

from ._cyclotomic import CyclotomicValue, _size
from .core import DomainError, Padic, PrecisionError, exact, p_power, valuation
from .harmonic import PadicTestFunction, _prime, _values, haar_volume


def kozyrev_wavelet(field, *, scale=0, center=0, frequency=1,
                    support_exponent=None, constancy_exponent=None, max_size=256):
    """chi(j*(x-center)/p**(scale+1)) on center+p**scale*Z_p, zero outside.

    frequency=j lies in 1..p-1. Values have amplitude one, mean zero and
    exact Haar norm squared p**(-scale). This avoids adjoining square roots
    for L2 normalization. Our scale is minus Kozyrev's gamma; centering the
    phase changes the published wavelet only by a constant unit phase.

    Optional lattice exponents embed the wavelet in a larger/finer grid.
    A p-adic centre needs absolute precision >=scale+1 to fix its phase.
    """
    p = _prime(field,center)
    if type(scale) is not int:
        raise TypeError("wavelet scale must be an integer")
    if type(frequency) is not int or not 1 <= frequency < p:
        raise ValueError("wavelet frequency must be an integer in 1..p-1")
    if isinstance(center,Padic):
        if center.precision < scale+1:
            raise PrecisionError("wavelet centre needs digits through scale+1")
        center = center.representative()
    center = exact(center)
    natural_a = min(scale,valuation(center,p))
    a = natural_a if support_exponent is None else support_exponent
    b = scale+1 if constancy_exponent is None else constancy_exponent
    if type(a) is not int or type(b) is not int:
        raise TypeError("lattice exponents must be integers")
    if a > natural_a or b < scale+1:
        raise DomainError("lattice must contain the wavelet support and resolve its phase")
    size = _size(p,b-a,max_size)
    step,unit = p_power(p,a),p_power(p,scale)
    roots = tuple(CyclotomicValue.root(p,1,frequency*j) for j in range(p))
    values = []
    for k in range(size):
        x = (step*k-center)/unit
        if valuation(x,p) < 0:
            values.append(0)
        else:
            residue = x.numerator*pow(x.denominator,-1,p) % p
            values.append(roots[residue])
    return PadicTestFunction(field,values,support_exponent=a,
                             constancy_exponent=b,max_size=max_size)


@dataclass(frozen=True, init=False)
class KozyrevBasis:
    """A finite orthogonal basis of a declared test-function lattice.

    The first mode is the container indicator. The others have labels
    (scale, translation_index, frequency), with centre=p**a*translation_index.
    Their exact squared norms are retained; the modes are not rescaled by
    possibly non-native square roots. There are exactly p**(b-a) modes.
    """
    template: PadicTestFunction
    labels: tuple
    modes: tuple[PadicTestFunction, ...]
    norms_squared: tuple[Fraction, ...]

    def __init__(self, field, *, support_exponent=0, constancy_exponent=1, max_size=256):
        p = _prime(field)
        a,b = support_exponent,constancy_exponent
        if type(a) is not int or type(b) is not int:
            raise TypeError("lattice exponents must be integers")
        if b < a:
            raise ValueError("constancy exponent must be >= support exponent")
        size = _size(p,b-a,max_size)
        template = PadicTestFunction(field,[0]*size,support_exponent=a,
                                    constancy_exponent=b,max_size=max_size)
        labels,modes,norms = [None],[template._replace_values([1]*size)],[haar_volume(p,a)]
        for scale in range(a,b):
            for r in range(p**(scale-a)):
                for j in range(1,p):
                    labels.append((scale,r,j))
                    modes.append(kozyrev_wavelet(field,scale=scale,center=p_power(p,a)*r,
                        frequency=j,support_exponent=a,constancy_exponent=b,max_size=max_size))
                    norms.append(haar_volume(p,scale))
        for name,value in (("template",template),("labels",tuple(labels)),
                           ("modes",tuple(modes)),("norms_squared",tuple(norms))):
            object.__setattr__(self,name,value)

    def decompose(self, function):
        """Exact inner-product coefficients divided by mode norms squared."""
        self.template._check(function)
        return tuple(function.inner_product(mode)/norm
                     for mode,norm in zip(self.modes,self.norms_squared))

    def _coefficients(self, coefficients):
        size = len(self.modes)
        coefficients = tuple(islice(iter(coefficients),size+1))
        if len(coefficients) != size:
            raise ValueError(f"basis requires exactly {size} coefficients")
        return _values(coefficients,self.template.field.p,size)

    def reconstruct(self, coefficients):
        """Exact synthesis from a complete coefficient vector."""
        coefficients = self._coefficients(coefficients)
        p = self.template.field.p
        values = [sum((c*mode.values[j] for c,mode in zip(coefficients,self.modes)),
                      CyclotomicValue(p)) for j in range(len(coefficients))]
        return self.template._replace_values(values)

    def energy(self, coefficients):
        """Exact weighted Parseval sum of coefficient norms times mode norms."""
        coefficients = self._coefficients(coefficients)
        return sum((c*c.conjugate()*norm for c,norm in zip(coefficients,self.norms_squared)),
                   CyclotomicValue(self.template.field.p))
