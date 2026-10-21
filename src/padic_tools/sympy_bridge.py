"""Optional symbolic function expressions with explicit p-adic evaluation.

Install ``padic-tools[symbolic]``. These nodes intentionally differ from
SymPy's complex zeta and theta, and do not imply analytic continuation between
the complex and p-adic domains.
"""
from __future__ import annotations

from fractions import Fraction
import sympy as sp
from sympy.core.function import ArgumentIndexError

from . import analytic, special
from .core import DomainError, Padic


class PadicExp(sp.Function):
    """Symbolic local p-adic exponential; domain checked on evaluation."""
    nargs = 1

    def fdiff(self, argindex=1):
        if argindex != 1:
            raise ArgumentIndexError(self, argindex)
        return self


class PadicLog(sp.Function):
    """Symbolic Iwasawa logarithm with log(p)=0."""
    nargs = 1

    def fdiff(self, argindex=1):
        if argindex != 1:
            raise ArgumentIndexError(self, argindex)
        return 1/self.args[0]


class PadicTheta3(sp.Function):
    """Symbolic multiplicative theta(q,z), with explicit two arguments."""
    nargs = 2


class PadicZeta(sp.Function):
    """Symbolic KL branch Z(s,p,r); p and r remain explicit metadata."""
    nargs = 3


def evaluate(expression, field, substitutions=None):
    """Evaluate rational operations and typed function nodes in a Q_p field.

    Substitution values can be exact rationals or existing Padic balls. Floats,
    unbound symbols, ordinary SymPy zeta, noninteger algebraic powers and
    implicit complex constants are rejected. All numeric domains are checked.
    """
    substitutions = {} if substitutions is None else substitutions
    expression = sp.sympify(expression)

    def walk(node):
        if node in substitutions:
            return field(substitutions[node])
        if isinstance(node, sp.Rational):
            return Fraction(int(node.p), int(node.q))
        if isinstance(node, sp.Symbol):
            raise ValueError(f"missing substitution for {node}")
        if isinstance(node, sp.Add):
            result = Fraction(0)
            for child in node.args:
                result += walk(child)
            return result
        if isinstance(node, sp.Mul):
            result = Fraction(1)
            for child in node.args:
                result *= walk(child)
            return result
        if isinstance(node, sp.Pow):
            if not isinstance(node.exp, sp.Integer):
                raise DomainError("noninteger powers require an explicit local branch")
            return walk(node.base)**int(node.exp)
        if isinstance(node, (PadicExp, PadicLog)):
            arg = field(walk(node.args[0]))
            return analytic.exp(arg) if isinstance(node, PadicExp) else analytic.log(arg)
        if isinstance(node, PadicTheta3):
            return special.theta3(field(walk(node.args[0])), walk(node.args[1]))
        if isinstance(node, PadicZeta):
            s, p, r = node.args
            if not isinstance(p, sp.Integer) or not isinstance(r, sp.Integer) or int(p) != field.p:
                raise DomainError("PadicZeta prime must match the numerical field and branch be integer")
            return special.zeta_p(walk(s), field, branch=int(r))
        raise DomainError(f"unsupported symbolic node: {node}; use the typed p-adic functions")

    result = walk(expression)
    return result if isinstance(result, Padic) else field(result)
