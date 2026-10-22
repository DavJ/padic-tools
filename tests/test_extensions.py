"""Independent polynomial/matrix oracles and ideal-ball containment checks."""
from fractions import Fraction as F
from itertools import product
import random

import pytest
import sympy as sp

from padic_tools import (
    DomainError, Polynomial, PrecisionError, Qp, QuadraticResidueElement,
    QuadraticResidueField, UnramifiedQuadraticElement, UnramifiedQuadraticField,
    valuation,
)
from padic_tools.core import p_power


CASES = [(2, [1, 1, 1]), (3, [1, 0, 1]), (5, [2, 0, 1]),
         (7, [1, 0, 1]), (3, [F(1, 4), 3, 1])]
t = sp.Symbol("t")


def expression(pair):
    return sum(sp.Rational(q.numerator, q.denominator)*t**i for i, q in enumerate(pair))


def oracle_pair(expression_, modulus):
    reduced = sp.Poly(sp.rem(expression_, modulus, t), t)
    return tuple(F(reduced.nth(i)) for i in (0, 1))


@pytest.mark.parametrize("p", [2, 3, 5, 7])
def test_constructor_certifies_irreducibility_by_independent_root_search(p):
    for b, c in product(range(p), repeat=2):
        reducible = any((j*j + b*j + c) % p == 0 for j in range(p))
        if reducible:
            with pytest.raises(DomainError):
                UnramifiedQuadraticField(Qp(p), [c, b, 1])
        else:
            L = UnramifiedQuadraticField(Qp(p), [c, b, 1])
            assert (L.degree, L.ramification_index, L.residue_degree) == (2, 1, 2)
            assert L.polynomial(L.generator).is_zero()
            assert L.uniformizer().valuation() == 1


def test_invalid_presentations_and_inexact_inputs_are_rejected():
    for f in ([], [1], [1, 1], [1, 0, 2], [1, 0, 0, 1], [F(1, 3), 0, 1]):
        with pytest.raises(DomainError):
            UnramifiedQuadraticField(Qp(3), f)
    with pytest.raises(DomainError):
        UnramifiedQuadraticField(Qp(2), [1, 0, 1])  # Q2(i) is ramified.
    with pytest.raises(TypeError):
        UnramifiedQuadraticField(3, [1, 0, 1])
    with pytest.raises(ValueError):
        QuadraticResidueField(9, [1, 0, 1])
    with pytest.raises(ValueError):
        UnramifiedQuadraticField(Qp(3), [1, 0, 1], name="")
    L = UnramifiedQuadraticField(Qp(3), [1, 0, 1])
    for q in (0.1, 1j, True):
        with pytest.raises(TypeError):
            L(1, q)
    with pytest.raises(TypeError):
        UnramifiedQuadraticElement(L, 1, 2)


@pytest.mark.parametrize("p,f", CASES)
def test_residue_field_exhaustive_arithmetic_and_frobenius(p, f):
    L = UnramifiedQuadraticField(Qp(p, 7), f)
    R = L.residue_field()
    elements = [R(a, b) for a, b in product(range(p), repeat=2)]
    assert len(set(elements)) == R.order == p*p
    assert Polynomial(f)(R.generator) == R(0)
    for x in elements:
        assert x**(p*p) == x
        assert x**p == x.frobenius() == x.conjugate()
        assert x + x.conjugate() == R(x.trace())
        assert x * x.conjugate() == R(x.norm())
        assert L.lift_residue(x).residue() == x
        if x:
            assert x*x.inverse() == R(1)
            assert x**(-2) == 1/(x*x)
        else:
            with pytest.raises(ZeroDivisionError):
                x.inverse()
        for y in elements:
            lx, ly = L.lift_residue(x), L.lift_residue(y)
            assert (lx+ly).residue() == x+y
            assert (lx*ly).residue() == x*y
            assert (x*y).frobenius() == x.frobenius()*y.frobenius()
    assert not bool(R(0)) and bool(R.generator)
    assert 2-R.generator == R(2)-R.generator


@pytest.mark.parametrize("p,f", CASES)
def test_local_arithmetic_against_sympy_quotient_field_and_multiplication_matrix(p, f):
    L = UnramifiedQuadraticField(Qp(p, 9), f)
    modulus = expression(tuple(map(F, f)))
    rng = random.Random(45258+p)
    for _ in range(22):
        a = tuple(F(rng.randint(-25, 25), rng.randint(1, 12)) for _ in range(2))
        b = tuple(F(rng.randint(1, 25), rng.randint(1, 12)) for _ in range(2))
        x, y = L(a), L(b)
        A, B = expression(a), expression(b)
        assert (x+y).contains(oracle_pair(A+B, modulus))
        assert (x-y).contains(oracle_pair(A-B, modulus))
        assert (x*y).contains(oracle_pair(A*B, modulus))
        assert (x/y).contains(oracle_pair(A*sp.invert(B, modulus, t), modulus))
        assert (x/b).contains(oracle_pair(A*sp.invert(B, modulus, t), modulus))
        assert (x+F(1, p)).contains(oracle_pair(A+sp.Rational(1, p), modulus))
        assert (F(1, p)*x).contains(tuple(q/p for q in a))
        assert (F(1, p)-x).contains((F(1, p)-a[0], -a[1]))
        if any(a):
            assert x.inverse().contains(oracle_pair(sp.invert(A, modulus, t), modulus))
            v = min(valuation(q, p) for q in a)
            assert x.valuation() == v
            assert x.relative_precision == 9
            assert x.absolute_value() == p_power(p, -v)
            assert x.norm().valuation() == 2*v
        # Determinant/trace of multiplication by A, independent of the
        # implementation's closed quadratic norm and trace formulas.
        matrix = sp.Matrix.hstack(*(sp.Matrix(oracle_pair(A*t**j, modulus)) for j in range(2)))
        assert x.norm().contains(F(matrix.det()))
        assert x.trace().contains(F(matrix.trace()))
        assert (x+x.conjugate()).is_congruent(x.trace())
        assert (x*x.conjugate()).is_congruent(x.norm())
        for n in (-2, 0, 1, 3):
            if n >= 0 or any(a):
                power = A**n if n >= 0 else sp.invert(A, modulus, t)**(-n)
                assert (x**n).contains(oracle_pair(power, modulus))


@pytest.mark.parametrize("p,f", CASES[:3])
@pytest.mark.parametrize("n", [-1, 1, 4])
def test_entire_input_balls_enclosed_after_perturbation(p, f, n):
    L = UnramifiedQuadraticField(Qp(p, 8), f)
    modulus = expression(tuple(map(F, f)))
    x = L(p_power(p, n-2), p_power(p, n-1), absprec=n)
    y = L(1+p_power(p, n-3), 2*p_power(p, n-3), absprec=n+1)
    h, k = p_power(p, n), p_power(p, n+1)
    for da, db in product((-2, 0, 1), repeat=2):
        a = tuple(q+d*h for q, d in zip(x.representative(), (da, db)))
        b = tuple(q+d*k for q, d in zip(y.representative(), (db, da)))
        A, B = expression(a), expression(b)
        assert (x+y).contains(oracle_pair(A+B, modulus))
        assert (x*y).contains(oracle_pair(A*B, modulus))
        assert x.inverse().contains(oracle_pair(sp.invert(A, modulus, t), modulus))
        assert (x/y).contains(oracle_pair(A*sp.invert(B, modulus, t), modulus))
        matrix = sp.Matrix.hstack(*(sp.Matrix(oracle_pair(A*t**j, modulus)) for j in range(2)))
        assert x.norm().contains(F(matrix.det()))
        assert x.trace().contains(F(matrix.trace()))


def test_precision_formulas_exact_scalars_and_zero_balls():
    L = UnramifiedQuadraticField(Qp(3, 12), [1, 0, 1], "i")
    x, y = L(9, 18, absprec=10), L(F(1, 3), 0, absprec=8)
    assert (x+y).precision == 8
    assert (x*y).precision == 9
    assert x.inverse().precision == 6
    assert (x*9).precision == 12
    assert (x/9).precision == 8
    assert (x/L.base(F(1, 3), absprec=8)).precision == 11
    assert (L.base(3, absprec=6)*x).precision == 8
    assert x.norm().precision == 12 and x.trace().precision == 10
    assert x.conjugate().precision == 10 and x.frobenius().conjugate() == x
    assert x+0 == x and x*1 == x and x/1 == x
    assert (x-x).is_zero() and (x-x).precision == 10
    assert (x*0).is_zero() and (x*0).precision == 10
    z = L.O(-2)
    assert z.v == -2 and z.relative_precision == 0
    assert z.contains(F(1, 9), F(1, 9)) and not z.contains(F(1, 27))
    assert z.norm().precision == -4 and z.norm().is_zero()
    assert (z*L.O(3)).precision == 1
    for operation in (z.valuation, z.absolute_value, z.inverse):
        with pytest.raises(PrecisionError):
            operation()
    with pytest.raises(TypeError):
        bool(x)
    with pytest.raises(ZeroDivisionError):
        x/0
    with pytest.raises(TypeError):
        x**True


def test_conversion_never_manufactures_digits_and_balls_have_structural_equality():
    K = Qp(3, 12)
    L = UnramifiedQuadraticField(K, [1, 0, 1])
    x = L(K(1, absprec=4), K(2, absprec=7))
    assert x.precision == x.a.precision == x.b.precision == 4
    same_field = UnramifiedQuadraticField(Qp(3, 25), [1, 0, 1], "i")
    assert same_field(x) == x and hash(same_field(x)) == hash(x)
    assert L(1, 2, absprec=6) != x
    assert len({x, same_field(x), L(1, 2, absprec=6)}) == 2
    assert x.is_congruent(L(1+3**4, 2, absprec=8))
    assert not x.is_congruent(L(2, 2))
    assert x.truncate(2).precision == 2
    for conversion in (lambda: L(x, absprec=5), lambda: x.truncate(5),
                       lambda: L(K(1, absprec=4), absprec=5),
                       lambda: x.is_congruent(L(1, 2), precision=5)):
        with pytest.raises(PrecisionError):
            conversion()
    for value in (Qp(5)(1), UnramifiedQuadraticField(K, [2, 1, 1])(1)):
        with pytest.raises(TypeError):
            L(value)
        with pytest.raises(TypeError):
            x*value
    with pytest.raises(TypeError):
        x.contains(x)
    for n in (True, 2.5):
        with pytest.raises(TypeError):
            L(1, absprec=n)
        with pytest.raises(TypeError):
            L.O(n)


def test_residue_domains_lifts_and_parent_checks():
    L = UnramifiedQuadraticField(Qp(3), [1, 0, 1])
    R = L.residue_field()
    assert L(4, 5).residue() == R(1, 2)
    assert L.O(1).residue() == R(0)
    assert L(3).residue() == R(0)
    with pytest.raises(DomainError):
        L(F(1, 3), 1).residue()
    with pytest.raises(PrecisionError):
        L(1, 1, absprec=0).residue()
    with pytest.raises(PrecisionError):
        L.lift_residue(R(1), absprec=0)
    with pytest.raises(DomainError):
        R(F(1, 3))
    with pytest.raises(TypeError):
        QuadraticResidueElement(R, True, 0)
    with pytest.raises(TypeError):
        R(QuadraticResidueField(5, [2, 0, 1])(1))
    with pytest.raises(TypeError):
        R.generator**F(1, 2)
    # Integer-coefficient lifting is a section of reduction, not an additive
    # map from a field of characteristic p into characteristic zero.
    assert not (L.lift_residue(R(2))+L.lift_residue(R(2))).is_congruent(L.lift_residue(R(1)), 2)
