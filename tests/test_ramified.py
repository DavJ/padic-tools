"""Eisenstein fields: independent quotient/matrix and whole-ball oracles."""
from fractions import Fraction as F
from itertools import product
import random

import pytest
import sympy as sp

from padic_tools import (
    DomainError, EisensteinQuadraticElement, EisensteinQuadraticField,
    Polynomial, PrecisionError, PrimeResidueElement, PrimeResidueField,
    Q2i, Qp, UnramifiedQuadraticField, valuation,
)
from padic_tools.core import p_power


CASES = [(2, [2, -2, 1], 2), (2, [-2, 0, 1], 3),
         (2, [2, 0, 1], 3), (3, [-3, 0, 1], 1),
         (5, [5, 5, 1], 1), (3, [F(3, 2), F(3, 5), 1], 1)]
t = sp.Symbol("t")


def expression(pair):
    return sum(sp.Rational(q.numerator, q.denominator)*t**i for i, q in enumerate(pair))


def reduced(expression_, modulus):
    poly = sp.Poly(sp.rem(expression_, modulus, t), t)
    return tuple(F(poly.nth(i)) for i in (0, 1))


def matrix(expression_, modulus):
    return sp.Matrix.hstack(*(sp.Matrix(reduced(expression_*t**j, modulus)) for j in (0, 1)))


def power(n, modulus):
    return sp.rem(t**n if n >= 0 else sp.invert(t, modulus, t)**(-n), modulus, t)


@pytest.mark.parametrize("p,f,d", CASES)
def test_eisenstein_data_valuation_normalization_and_different(p, f, d):
    L = EisensteinQuadraticField(Qp(p, 9), f)
    assert (L.degree, L.ramification_index, L.residue_degree) == (2, 2, 1)
    assert L.precision == 18 and L.different_exponent == d
    assert valuation(F(f[1])**2-4*F(f[0]), p) == d
    pi = L.uniformizer()
    assert pi == L.generator and pi.valuation() == 1
    assert pi.base_valuation() == F(1, 2)
    assert L(p).valuation() == 2
    assert pi.absolute_value_squared() == F(1, p)
    assert L(p).absolute_value_squared() == F(1, p*p)
    assert Polynomial(f)(pi).is_zero()
    assert (pi*pi.inverse()).is_congruent(1)
    assert L.residue_field().order == p
    assert pi.residue() == L.residue_field()(0)
    assert not hasattr(pi, "frobenius")  # A ramified lift is not selected.


def test_validation_rejects_unverified_and_non_eisenstein_bases():
    for f in ([], [2], [2, 1], [2, 0, 2], [2, 0, 0, 1],
              [1, 0, 1], [4, 0, 1], [2, 1, 1], [F(1, 2), 0, 1]):
        with pytest.raises(DomainError):
            EisensteinQuadraticField(Qp(2), f)
    with pytest.raises(TypeError):
        EisensteinQuadraticField(2, [2, -2, 1])
    with pytest.raises(TypeError):
        EisensteinQuadraticElement(Q2i(), 1, 2)
    with pytest.raises(ValueError):
        Q2i(0)
    with pytest.raises(ValueError):
        EisensteinQuadraticField(Qp(2), [2, -2, 1], name="1pi")
    for scalar in (True, 0.1, 1j):
        with pytest.raises(TypeError):
            Q2i()(scalar)


def test_q2i_uses_uniformizer_basis_and_detects_gaussian_cancellation():
    L = Q2i(12)
    pi, i = L.generator, L.generator-1
    assert pi.precision == 25 and pi.relative_precision == 24
    assert L.polynomial == Polynomial([2, -2, 1])
    assert (i*i).is_congruent(-1)
    assert (1+i).valuation() == (1-i).valuation() == 1
    assert (pi*pi).is_congruent(2*i)
    assert pi.norm().contains(2) and pi.trace().contains(2)
    assert pi.conjugate().is_congruent(2-pi)
    assert i.conjugate().is_congruent(-i)
    assert L(F(1, 2)).valuation() == -2
    assert (1/pi).valuation() == -1
    # In the Gaussian basis, 1+i has two unit coefficients but valuation one.
    assert L(0, 1).contains(0, 1)
    assert not L(0, 1).contains(1, 1)


@pytest.mark.parametrize("p,f,d", CASES)
def test_exact_arithmetic_against_sympy_quotients_and_matrices(p, f, d):
    L = EisensteinQuadraticField(Qp(p, 8), f)
    modulus = expression(tuple(map(F, f)))
    rng = random.Random(8125+p+d)
    for _ in range(18):
        a = tuple(F(rng.randint(-20, 20), rng.randint(1, 14)) for _ in (0, 1))
        b = tuple(F(rng.randint(1, 20), rng.randint(1, 14)) for _ in (0, 1))
        A, B = expression(a), expression(b)
        x, y = L(a), L(b)
        assert (x+y).contains(reduced(A+B, modulus))
        assert (x-y).contains(reduced(A-B, modulus))
        assert (x*y).contains(reduced(A*B, modulus))
        assert (x/y).contains(reduced(A*sp.invert(B, modulus, t), modulus))
        assert (x/b).contains(reduced(A*sp.invert(B, modulus, t), modulus))
        assert (F(1, p)-x).contains((F(1, p)-a[0], -a[1]))
        assert (Qp(p)(1)+x).contains((a[0]+1, a[1]))
        M = matrix(A, modulus)
        assert x.norm().contains(F(M.det()))
        assert x.trace().contains(F(M.trace()))
        assert (x+x.conjugate()).is_congruent(x.trace())
        assert (x*x.conjugate()).is_congruent(x.norm())
        if any(a):
            # The residue degree is one, hence v_p(det(m_x))=v_pi(x).
            v = valuation(F(M.det()), p)
            assert x.valuation() == v and x.base_valuation() == F(v, 2)
            assert x.relative_precision == 16
            assert x.norm().valuation() == v
            assert x.inverse().contains(reduced(sp.invert(A, modulus, t), modulus))
        for n in (-2, 0, 1, 3):
            if n >= 0 or any(a):
                P = A**n if n >= 0 else sp.invert(A, modulus, t)**(-n)
                assert (x**n).contains(reduced(P, modulus))


@pytest.mark.parametrize("p,f,d", [CASES[0], CASES[1], CASES[3]])
@pytest.mark.parametrize("n", [-3, -2, 1, 2, 5, 6])
def test_perturbing_inputs_by_entire_uniformizer_ideals(p, f, d, n):
    L = EisensteinQuadraticField(Qp(p, 8), f)
    modulus = expression(tuple(map(F, f)))
    A0, B0 = power(n-2, modulus)*(1+t), power(n-3, modulus)*(1+t)
    x = L(reduced(A0, modulus), absprec=n)
    y = L(reduced(B0, modulus), absprec=n+2)
    assert [c.precision for c in x.coefficients] == [(n+1)//2, n//2]
    for da, db in product((-1, 0, 1), repeat=2):
        A = A0+power(n, modulus)*(da+db*t)
        B = B0+power(n+2, modulus)*(db+da*t)
        assert x.contains(reduced(A, modulus))
        assert y.contains(reduced(B, modulus))
        assert (x+y).contains(reduced(A+B, modulus))
        assert (x*y).contains(reduced(A*B, modulus))
        assert x.inverse().contains(reduced(sp.invert(A, modulus, t), modulus))
        assert (x/y).contains(reduced(A*sp.invert(B, modulus, t), modulus))
        M = matrix(A, modulus)
        assert x.norm().contains(F(M.det()))
        assert x.trace().contains(F(M.trace()))


@pytest.mark.parametrize("n", range(-5, 8))
def test_coordinate_staggering_and_zero_ball_bounds(n):
    L = Q2i(8)
    z = L.O(n)
    assert z.is_zero() and z.v == n and z.relative_precision == 0
    assert z.a.precision == (n+1)//2 and z.b.precision == n//2
    assert z.norm().precision == n and z.norm().is_zero()
    assert z.trace().precision == (n+2)//2 and z.trace().is_zero()
    assert z.contains(p_power(2, (n+1)//2), 0)
    assert z.contains(0, p_power(2, n//2))
    assert not z.contains(p_power(2, (n+1)//2-1), 0)
    assert not z.contains(0, p_power(2, n//2-1))
    for op in (z.valuation, z.base_valuation, z.absolute_value_squared, z.inverse):
        with pytest.raises(PrecisionError):
            op()
    assert (z*L.O(2)).precision == n+2
    if n >= 1:
        assert z.residue() == L.residue_field()(0)
    else:
        with pytest.raises((DomainError, PrecisionError)):
            z.residue()


def test_precision_arithmetic_and_coordinate_input_caps():
    L = Q2i(8)
    x, y = L(0, 2, absprec=11), L(F(1, 2), absprec=7)
    assert x.v == 3 and y.v == -2
    assert (x+y).precision == 7
    assert (x*y).precision == 9
    assert x.inverse().precision == 5
    assert (x*2).precision == 13 and (x/2).precision == 9
    assert x.trace().precision == 6 and x.norm().precision == 8
    assert (x-x).precision == 11 and (x-x).is_zero()
    assert x.conjugate().conjugate() == x
    assert x+0 == x and x*1 == x and x/1 == x
    assert (x*0).is_zero() and (x*0).precision == 11
    assert L(Qp(2)(1, absprec=4)).precision == 8
    mixed = L(Qp(2)(1, absprec=4), Qp(2)(2, absprec=3))
    assert mixed.precision == 7
    assert [c.precision for c in mixed.coefficients] == [4, 3]
    assert (Qp(2)(2, absprec=4)*x).precision == 11
    for op in (lambda: L(mixed, absprec=8), lambda: mixed.truncate(8),
               lambda: L(0, Qp(2)(1, absprec=3), absprec=8),
               lambda: mixed.is_congruent(L(1, 2), precision=8)):
        with pytest.raises(PrecisionError):
            op()
    assert mixed.truncate(5).precision == 5
    for op in (lambda: L.O(True), lambda: L(1, absprec=1.5)):
        with pytest.raises(TypeError):
            op()
    with pytest.raises(TypeError):
        bool(x)
    with pytest.raises(TypeError):
        x**F(1, 2)
    with pytest.raises(ZeroDivisionError):
        x/0


def test_exact_small_quotient_sizes_and_structural_equality():
    L = Q2i(8)
    for n in range(1, 7):
        residues = {L(a, b, absprec=n) for a in range(2**((n+1)//2))
                    for b in range(2**(n//2))}
        assert len(residues) == 2**n
        for x in residues:
            assert x.conjugate().conjugate() == x
            if not x.is_zero():
                assert (x*x.inverse()).is_congruent(1)
    x = L(1, 2, absprec=5)
    same = EisensteinQuadraticField(Qp(2, 20), [2, -2, 1], name="u")(x)
    assert same == x and hash(same) == hash(x)
    assert x != x.truncate(4) and x.is_congruent(x.truncate(4))
    assert len({x, same, x.truncate(4)}) == 2
    for foreign in (UnramifiedQuadraticField(Qp(2), [1, 1, 1])(1),
                    EisensteinQuadraticField(Qp(2), [-2, 0, 1])(1), Qp(3)(1)):
        with pytest.raises(TypeError):
            L(foreign)
        with pytest.raises(TypeError):
            x*foreign
    with pytest.raises(TypeError):
        x.contains(x)


@pytest.mark.parametrize("p", [2, 3, 5, 7])
def test_prime_residue_field_exact_arithmetic_and_reduction(p):
    R = PrimeResidueField(p)
    L = EisensteinQuadraticField(Qp(p, 8), [p, 0, 1])
    assert R.order == p and len({R(a) for a in range(p)}) == p
    for a in range(p):
        x = R(a)
        assert x**p == x.frobenius() == x
        assert x.trace() == x.norm() == a
        assert L.lift_residue(x).residue() == x
        assert L(a, 1, absprec=1).residue() == x
        assert L(a, 1).conjugate().residue() == x
        if a:
            assert x*x.inverse() == R(1)
            assert x**(-2) == 1/(x*x)
        else:
            with pytest.raises(ZeroDivisionError):
                x.inverse()
        for b in range(p):
            y = R(b)
            assert (L(a, 1)+L(b, 2)).residue() == x+y
            assert (L(a, 1)*L(b, 2)).residue() == x*y
            assert (x+y).value == (a+b) % p
            assert (x*y).value == a*b % p
        assert 2-x == R(2)-x
    assert R(-1).value == p-1
    assert bool(R(1)) and not bool(R(0))


def test_residue_domain_errors_and_lifts():
    L = Q2i()
    R = L.residue_field()
    with pytest.raises(DomainError):
        (1/L.generator).residue()
    with pytest.raises(PrecisionError):
        L(1, absprec=0).residue()
    with pytest.raises(PrecisionError):
        L.lift_residue(R(1), absprec=0)
    with pytest.raises(DomainError):
        R(F(1, 2))
    with pytest.raises(TypeError):
        R(PrimeResidueField(3)(1))
    with pytest.raises(TypeError):
        PrimeResidueElement(R, True)
    with pytest.raises(TypeError):
        R(1)**True
    with pytest.raises(ValueError):
        PrimeResidueField(9)
    assert not (L.lift_residue(R(1))+L.lift_residue(R(1))).is_congruent(L.lift_residue(R(0)), 3)
