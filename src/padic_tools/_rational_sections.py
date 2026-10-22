"""Exact residue-class sections of rational generating functions over Q.

Uses decimated linear recurrences and Cayley--Hamilton, without adjoining
roots of unity. This is an internal rational polynomial algorithm.
"""
from fractions import Fraction

from .core import PrecisionError
from .symbolic import Polynomial


def _identity(n):
    return [[Fraction(i == j) for j in range(n)] for i in range(n)]


def _multiply(a, b):
    columns = list(zip(*b))
    return [[sum((x*y for x,y in zip(row,col)),Fraction(0))
             for col in columns] for row in a]


def _power(a, n):
    result = _identity(len(a))
    while n:
        if n & 1:
            result = _multiply(result,a)
        n >>= 1
        if n:
            a = _multiply(a,a)
    return result


def _vector(a, v):
    return [sum((x*y for x,y in zip(row,v)),Fraction(0)) for row in a]


def _characteristic(a):
    # Faddeev--LeVerrier: coefficients of det(I-u*A), constant first.
    n = len(a)
    b, cs = _identity(n),[Fraction(1)]
    for k in range(1,n+1):
        b = _multiply(a,b)
        c = -sum(b[i][i] for i in range(n))/k
        cs.append(c)
        for i in range(n):
            b[i][i] += c
    return Polynomial(cs)


def _inflate(poly, step):
    if poly.degree < 0:
        return poly
    cs = [Fraction(0)]*(poly.degree*step+1)
    for n,c in enumerate(poly.coefficients):
        cs[n*step] = c
    return Polynomial(cs)


def zero_section(numerator, denominator, period, max_degree):
    """Keep powers z**n with n divisible by period in a rational Laurent series.

    The expansion is at z=0, allowing a finite pole. Return exact rational
    numerator/denominator data in z, with no series truncation.
    """
    m = next(n for n,c in enumerate(denominator.coefficients) if c)
    regular = Polynomial(denominator.coefficients[m:])
    polynomial,proper = divmod(numerator,regular)
    d,r,offset = regular.degree,m % period,(m//period)*period
    # Size limits concern this exact native algorithm, not mathematical validity.
    degree = max(offset, max(0,polynomial.degree-r)) + period*d
    if degree > max_degree:
        raise PrecisionError(f"unit restriction may require degree {degree}; max_degree={max_degree}")
    if d > 64:
        raise PrecisionError("native unit restriction supports recurrence order at most 64")
    selected = Polynomial()
    if polynomial.degree >= r:
        cs = [Fraction(0)]*(polynomial.degree-r+1)
        for n in range(r,polynomial.degree+1,period):
            cs[n-r] = polynomial[n]
        selected = Polynomial(cs)
    if d and proper.degree >= 0:
        q = regular*(1/regular[0])
        p = proper*(1/regular[0])
        initial = []
        for n in range(d):
            initial.append(p[n]-sum((q[j]*initial[n-j] for j in range(1,n+1)),Fraction(0)))
        companion = [[Fraction(0)]*d for _ in range(d)]
        for i in range(d-1):
            companion[i][i+1] = Fraction(1)
        companion[-1] = [-q[d-i] for i in range(d)]
        stride = _power(companion,period)
        state = _vector(_power(companion,r),initial)
        values = []
        for _ in range(d):
            values.append(state[0])
            state = _vector(stride,state)
        decimated_denominator = _characteristic(stride)
        decimated_numerator = Polynomial([
            sum((decimated_denominator[j]*values[n-j] for j in range(n+1)),Fraction(0))
            for n in range(d)])
        bottom = _inflate(decimated_denominator,period)
        top = _inflate(decimated_numerator,period)
    else:
        top,bottom = Polynomial(),Polynomial([1])
    top = selected*bottom+top
    if offset:
        bottom = Polynomial([0]*offset+list(bottom.coefficients))
    return top,bottom
