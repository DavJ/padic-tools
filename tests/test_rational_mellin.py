from fractions import Fraction as F

import pytest

from padic_tools import (
    AtomicMeasure, DomainError, PrecisionError, Qp,
    RationalAmiceMeasure, power, teichmuller, zeta_p,
)


@pytest.mark.parametrize("p", [2,3,5,7])
def test_exact_unit_restriction_matches_mixed_integer_atomic_support(p):
    K = Qp(p,12)
    atom = AtomicMeasure(K,[(-2*p,2),(-p+1,F(1,p)),(-1,3),(0,7),
                           (1,-2),(p,4),(p+1,3),(2*p,1)])
    source = RationalAmiceMeasure.from_atomic(atom)
    units = source.restrict_units()
    expected = RationalAmiceMeasure.from_atomic(atom.restrict_units())
    assert units == expected
    assert units.restrict_units() is units
    for n in range(10):
        assert units.moment(n) == atom.restrict_units().moment(n)
    assert units.amice(p).is_congruent(atom.restrict_units().amice(p))
    with pytest.raises(DomainError): source.mellin(1)


@pytest.mark.parametrize("p", [2,3,5])
@pytest.mark.parametrize("numerator,denominator", [([1,2],[3,1,1]),([1],[1,3,3,1])])
def test_rational_restriction_against_independent_cyclotomic_quotient(p,numerator,denominator):
    sp = pytest.importorskip('sympy')
    K = Qp(p,12)
    if denominator[0] % p == 0:
        denominator = [denominator[0]+1,*denominator[1:]]
    mu = RationalAmiceMeasure(K,numerator,denominator)
    restricted = mu.restrict_units()
    z,w = sp.symbols('z w')
    domain = sp.QQ.frac_field(z)
    # Independent method: inversion in Q(z)[w]/(w**p-1). Its constant
    # coefficient is the root-of-unity average, without approximating roots.
    top = mu.numerator.to_sympy('t').subs(sp.Symbol('t'),z*w-1)
    bottom = mu.denominator.to_sympy('t').subs(sp.Symbol('t'),z*w-1)
    modulus = sp.Poly(w**p-1,w,domain=domain)
    inverse = sp.invert(sp.Poly(bottom,w,domain=domain),modulus)
    projected = (sp.Poly(top,w,domain=domain)*inverse).rem(modulus).nth(0)
    original = mu.numerator.to_sympy('z')/mu.denominator.to_sympy('z')
    actual = restricted.numerator.to_sympy('z')/restricted.denominator.to_sympy('z')
    assert sp.cancel(actual-original+projected.subs(z,z+1)) == 0


@pytest.mark.parametrize("p", [3,5,7,137,139])
def test_restricted_fermionic_moments_have_exact_euler_factor(p):
    K = Qp(p,8)
    source = RationalAmiceMeasure(K,1,[2,1])
    units = source.restrict_units()
    for n in range(9):
        assert units.moment(n) == (1-p**n)*source.moment(n)
    # Restriction already proved unit support: idempotence needs no larger
    # recurrence, even at the UBT primes 137 and 139.
    assert units.restrict_units(max_degree=1) is units


@pytest.mark.parametrize("p", [2,3,5,7])
def test_mellin_matches_atomic_noninteger_and_integer_character_values(p):
    K = Qp(p,12)
    atomic = AtomicMeasure(K,[(-1,F(1,p)),(1,2),(p+1,-3),(2*p+1,4)])
    mu = RationalAmiceMeasure.from_atomic(atomic)
    period = 2 if p == 2 else p-1
    for s in [F(1,3 if p == 2 else 2),F(-1),F(0),F(7)]:
        for branch in range(period):
            value = mu.mellin(s,branch=branch)
            assert value.precision == 12
            assert value.is_congruent(atomic.mellin(s,branch=branch),precision=12)
    for n in range(9):
        assert mu.mellin(n,branch=n % period).contains(atomic.moment(n))


@pytest.mark.parametrize("p", [3,5,7])
def test_nonatomic_mellin_matches_independent_kl_zeta_algorithm(p):
    K = Qp(p,12)
    mu = RationalAmiceMeasure(K,1,[2,1]).restrict_units()
    s = F(1,2)
    principal2 = K(2)/teichmuller(K(2))
    expected = (power(principal2,s+1)-1)*zeta_p(-s,K,branch=0)
    value = mu.mellin(s,branch=p-2)
    assert value.precision == 12
    assert value.is_congruent(expected)
    for n in range(8):
        assert mu.mellin(n,branch=n % (p-1)).contains(mu.moment(n))
    # The removable KL-expression singularity at s=-1 is evaluated directly.
    assert mu.mellin(-1,branch=p-2).precision == 12
    # Symmetry forces even character branches to vanish for this measure.
    assert mu.mellin(s,branch=0).is_zero()


def test_dyadic_nonatomic_mellin_interpolates_exact_restricted_moments():
    mu = RationalAmiceMeasure(Qp(2,12),1,[3,1]).restrict_units()
    for n in range(16):
        value = mu.mellin(n,branch=n % 2)
        assert value.precision == 12
        assert value.contains(mu.moment(n))
    for branch in [0,1]:
        assert mu.mellin(F(1,3),branch=branch).precision == 12


@pytest.mark.parametrize("p", [2,3,5])
def test_exponent_precision_gain_and_absolute_requested_precision(p):
    K = Qp(p,10)
    atom = AtomicMeasure(K,[(1,F(1,p*p)),(p+1,2)])
    mu = RationalAmiceMeasure.from_atomic(atom)
    q = F(1,3 if p == 2 else 2)
    s = K(q,absprec=3)
    gain = 2 if p == 2 else 1
    available = 3+gain-2
    value = mu.mellin(s,branch=1)
    assert value.precision == available
    assert value.is_congruent(atom.mellin(q,branch=1),precision=available)
    with pytest.raises(PrecisionError): mu.mellin(s,branch=1,precision=available+1)
    # Exact source data may be evaluated beyond the default field precision.
    value = mu.mellin(q,branch=1,precision=16)
    assert value.precision == 16
    reference = AtomicMeasure(Qp(p,25),atom.atoms).mellin(q,branch=1)
    assert value.is_congruent(reference,precision=16)
    assert mu.mellin(K.O(0),branch=1).precision == gain-2


def test_unit_support_proofs_follow_addition_scaling_but_not_additive_convolution():
    K = Qp(3,10)
    a = RationalAmiceMeasure.from_atomic(AtomicMeasure(K,[(1,1)]))
    b = RationalAmiceMeasure.from_atomic(AtomicMeasure(K,[(2,1)]))
    assert (a+b).scale(3).mellin(0).contains(6)
    with pytest.raises(DomainError): a.convolve(b).mellin(0)
    assert a.convolve(b).restrict_units().mellin(0).is_zero()
    assert RationalAmiceMeasure(K,0).mellin(0).is_zero()
    # The proof flag is metadata, not part of mathematical measure equality.
    assert a == RationalAmiceMeasure(K,[1,1])


def test_resource_limits_domains_and_type_validation_are_explicit():
    K = Qp(5,12)
    source = RationalAmiceMeasure(K,1,[2,1])
    with pytest.raises(PrecisionError): source.restrict_units(max_degree=1)
    for invalid in [0,-1,True]:
        with pytest.raises(ValueError): source.restrict_units(max_degree=invalid)
    units = source.restrict_units()
    with pytest.raises(PrecisionError): units.mellin(F(1,2),max_terms=2)
    with pytest.raises(DomainError): units.mellin(F(1,5))
    with pytest.raises(TypeError): units.mellin(Qp(3)(1))
    with pytest.raises(TypeError): units.mellin(1,branch=True)
    with pytest.raises(TypeError): units.mellin(1,precision=True)
    with pytest.raises(ValueError): units.mellin(1,max_terms=True)
    # A tiny measure is known to vanish at the requested precision.
    tiny = RationalAmiceMeasure.from_atomic(AtomicMeasure(K,[(1,5**30)]))
    assert tiny.mellin(1,precision=12).is_zero()
    assert units.mellin(F(1,2),precision=0).precision == 0
