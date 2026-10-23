"""Dense exact and local-field matrices with certified pivot precision.

Local matrices are Cartesian products of scalar balls. Elimination uses
complete pivoting by the smallest valuation, including lower bounds from
zero-containing balls. Such balls are never silently treated as exact zero.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from .algebra import Gaussian, _exact_zero
from .bridges import GaussianBridge
from .core import DomainError, PadicField, PrecisionError, exact
from .extensions import UnramifiedQuadraticField
from .ramified import EisensteinQuadraticField


_LOCAL_FIELDS = (PadicField, UnramifiedQuadraticField, EisensteinQuadraticField)
_MAX_DIMENSION = 64


def _check_base(base):
    if base is not exact and base != Gaussian and not isinstance(base, _LOCAL_FIELDS):
        raise TypeError("base must be exact, Gaussian, or a certified local field")


def _convert(base, value, absprec=None):
    if isinstance(base, _LOCAL_FIELDS):
        return base(value, absprec=absprec)
    if absprec is not None:
        raise TypeError("absprec applies only to local-field matrices")
    return exact(value) if base is exact else GaussianBridge(value).value


def _multiply(a, b):
    # Internal structural zeros (e.g. an exact RHS) have no uncertainty.
    # Scalar zero balls do have uncertainty, so _exact_zero excludes them.
    return Fraction(0) if _exact_zero(a) or _exact_zero(b) else a*b


def _dot(a, b):
    return sum((_multiply(x, y) for x, y in zip(a, b)), Fraction(0))


def _zero_possible(value):
    return value.is_zero() if hasattr(value, "is_zero") else _exact_zero(value)


@dataclass(frozen=True)
class LinearSolveResult:
    """A solution enclosure and the certified elimination diagnostics.

    Permutations list original row/column indices in elimination order.
    For local fields, ``pivot_valuations`` are nondecreasing elementary
    divisor valuations in uniformizer units. For exact bases they are empty,
    and condition exponents are None. ``precisions`` mirrors the solution
    shape, using None for exact entries.
    """
    solution: object
    row_permutation: tuple
    column_permutation: tuple
    pivot_valuations: tuple
    precisions: tuple
    ramification_index: int | None

    @property
    def condition_exponent(self):
        """s_max-s_min in uniformizer units, or None for an exact base."""
        if self.ramification_index is None:
            return None
        return self.pivot_valuations[-1]-self.pivot_valuations[0]

    @property
    def base_condition_exponent(self):
        """The exact exponent of p in the matrix sup-norm condition number."""
        if self.ramification_index is None:
            return None
        return Fraction(self.condition_exponent, self.ramification_index)


@dataclass(frozen=True, init=False)
class Matrix:
    """Immutable nonempty dense matrix over Q, Q[i], or a local field.

    ``Matrix(exact, rows)`` and ``Matrix(Gaussian, rows)`` keep exact sources.
    With a local field, inputs become scalar balls. ``absprec`` optionally
    sets every entry's absolute uniformizer precision; existing balls cannot
    gain digits. Each dimension is limited to 64 in this native backend.

    Matrix products use ``@``; ``*`` is scalar multiplication. Bases must
    agree for arithmetic. Use ``map_coefficients`` for an explicit base change.
    """
    base: object
    rows: tuple

    def __init__(self, base, rows, *, absprec=None):
        _check_base(base)
        if absprec is not None and type(absprec) is not int:
            raise TypeError("absprec must be an integer")
        source = tuple(tuple(row) for row in rows)
        if not source or not source[0]:
            raise ValueError("matrices must have at least one row and column")
        ncols = len(source[0])
        if any(len(row) != ncols for row in source):
            raise ValueError("matrix rows must have equal length")
        if max(len(source), ncols) > _MAX_DIMENSION:
            raise DomainError("native dense matrices support dimensions <= 64")
        values = tuple(tuple(_convert(base, c, absprec) for c in row) for row in source)
        object.__setattr__(self, "base", base)
        object.__setattr__(self, "rows", values)

    @classmethod
    def identity(cls, base, size, *, absprec=None):
        """Identity entries, rounded to balls when the base is a local field."""
        if type(size) is not int or not 1 <= size <= _MAX_DIMENSION:
            raise ValueError("identity size must be an integer from 1 to 64")
        return cls(base, ((int(i == j) for j in range(size)) for i in range(size)), absprec=absprec)

    @property
    def nrows(self):
        """Number of rows."""
        return len(self.rows)

    @property
    def ncols(self):
        """Number of columns."""
        return len(self.rows[0])

    @property
    def shape(self):
        """The pair (nrows, ncols)."""
        return self.nrows, self.ncols

    @property
    def precisions(self):
        """Absolute uniformizer precision of each entry; None for exact entries."""
        return tuple(tuple(getattr(c, "precision", None) for c in row) for row in self.rows)

    def __getitem__(self, index):
        if isinstance(index, tuple):
            i, j = index
            return self.rows[i][j]
        return self.rows[index]

    def __iter__(self):
        return iter(self.rows)

    def _compatible(self, other):
        if not isinstance(other, Matrix) or self.base != other.base:
            raise TypeError("matrix arithmetic requires the same base; use map_coefficients")

    def _same_shape(self, other):
        self._compatible(other)
        if self.shape != other.shape:
            raise ValueError("matrix shapes must agree")

    def _square(self):
        if self.nrows != self.ncols:
            raise ValueError("operation requires a square matrix")

    def __add__(self, other):
        self._same_shape(other)
        return Matrix(self.base, ((a+b for a, b in zip(x, y)) for x, y in zip(self.rows, other.rows)))

    def __neg__(self):
        return Matrix(self.base, ((-c for c in row) for row in self.rows))

    def __sub__(self, other):
        return self+-other

    def _scalar(self, value):
        if isinstance(self.base, _LOCAL_FIELDS):
            try:
                # Exact rational scaling must not introduce default rounding.
                return exact(value)
            except TypeError:
                return self.base(value)
        return _convert(self.base, value)

    def __mul__(self, scalar):
        if isinstance(scalar, Matrix):
            return NotImplemented
        c = self._scalar(scalar)
        return Matrix(self.base, ((x*c for x in row) for row in self.rows))

    __rmul__ = __mul__

    def __matmul__(self, other):
        self._compatible(other)
        if self.ncols != other.nrows:
            raise ValueError("inner matrix dimensions must agree")
        columns = tuple(zip(*other.rows))
        return Matrix(self.base, ((_dot(row, column) for column in columns) for row in self.rows))

    def __pow__(self, exponent):
        self._square()
        if type(exponent) is not int:
            raise TypeError("matrix exponent must be an integer")
        if exponent < 0:
            return self.inverse()**(-exponent)
        if exponent == 0:
            return Matrix.identity(self.base, self.nrows)
        result, factor = None, self
        while exponent:
            if exponent & 1:
                result = factor if result is None else result@factor
            exponent >>= 1
            if exponent:
                factor = factor@factor
        return result

    def transpose(self):
        """Transpose, without scalar conjugation."""
        return Matrix(self.base, zip(*self.rows))

    @property
    def T(self):
        """The transpose."""
        return self.transpose()

    def trace(self):
        """Sum of the diagonal entries of a square matrix."""
        self._square()
        return sum((self.rows[i][i] for i in range(self.nrows)), Fraction(0))

    def determinant(self):
        """Division-free Berkowitz determinant, including unresolved zero balls.

        No pivot test is used. A local zero-containing result does not certify
        exact singularity. The native arithmetic cost is O(n**4).
        """
        self._square()
        coefficients = [Fraction(1)]  # Descending coefficients of det(tI-B).
        for start in range(self.nrows-1, -1, -1):
            m = self.nrows-start
            q = [Fraction(1), -self.rows[start][start]]
            row = self.rows[start][start+1:]
            vector = [self.rows[i][start] for i in range(start+1, self.nrows)]
            for k in range(m-1):
                q.append(-_dot(row, vector))
                if k < m-2:
                    vector = [_dot(self.rows[i][start+1:], vector)
                              for i in range(start+1, self.nrows)]
            coefficients = [sum((_multiply(q[i-j], coefficients[j])
                                 for j in range(min(i+1, len(coefficients)))), Fraction(0))
                            for i in range(m+1)]
        return coefficients[-1]*(-1)**self.nrows

    def map_coefficients(self, target_base, mapper=None, *, absprec=None):
        """Apply an explicit scalar map, e.g. an exact Gaussian local place."""
        _check_base(target_base)
        if mapper is None:
            return Matrix(target_base, self.rows, absprec=absprec)
        return Matrix(target_base, ((mapper(c) for c in row) for row in self.rows), absprec=absprec)

    def truncate(self, precision):
        """Truncate every local entry to one absolute precision; never add digits."""
        if not isinstance(self.base, _LOCAL_FIELDS):
            raise TypeError("exact matrices do not have finite precision")
        return Matrix(self.base, self.rows, absprec=precision)

    def is_congruent(self, other, precision=None):
        """Entrywise known-digit comparison; bases and shapes must agree."""
        self._same_shape(other)
        if precision is not None and type(precision) is not int:
            raise TypeError("precision must be an integer")
        if not isinstance(self.base, _LOCAL_FIELDS):
            if precision is not None:
                raise TypeError("exact matrices have no congruence precision")
            return self == other
        return all(a.is_congruent(b, precision) for x, y in zip(self.rows, other.rows) for a, b in zip(x, y))

    def contains(self, rows):
        """Test entrywise membership of an exact matrix or exact scalar rows."""
        source = rows.rows if isinstance(rows, Matrix) else tuple(tuple(row) for row in rows)
        if len(source) != self.nrows or any(len(row) != self.ncols for row in source):
            raise ValueError("membership requires matching shapes")
        if not isinstance(self.base, _LOCAL_FIELDS):
            return self == Matrix(self.base, source)
        return all(a.contains(b) for x, y in zip(self.rows, source) for a, b in zip(x, y))

    def solve(self, rhs):
        """Solve A*x=rhs, returning a tuple for a vector or Matrix for many RHS.

        Local results enclose the solutions for all represented input values.
        A missing or unresolved minimum-valuation pivot raises PrecisionError;
        exact singular systems raise ZeroDivisionError. Rectangular or singular
        solution families are outside this solver's scope.
        """
        return self.solve_with_info(rhs).solution

    def solve_with_info(self, rhs):
        """Solve with certified row/column permutations and condition exponents."""
        self._square()
        if isinstance(rhs, Matrix):
            self._compatible(rhs)
            if rhs.nrows != self.nrows:
                raise ValueError("right-hand side row count must match")
            b = [list(row) for row in rhs.rows]
            vector_rhs = False
        else:
            values = tuple(rhs)
            if len(values) != self.nrows:
                raise ValueError("right-hand side length must match")
            b = [[self._scalar(c)] for c in values]
            vector_rhs = True
        return self._eliminate(b, vector_rhs)

    def inverse(self):
        """Certified inverse, solving against an internally exact identity RHS."""
        self._square()
        b = [[Fraction(int(i == j)) for j in range(self.nrows)] for i in range(self.nrows)]
        return self._eliminate(b, False).solution

    def _eliminate(self, b, vector_rhs):
        local = isinstance(self.base, _LOCAL_FIELDS)
        n, width = self.nrows, len(b[0])
        a = [list(row) for row in self.rows]
        row_order, column_order, pivots = list(range(n)), list(range(n)), []
        for k in range(n):
            candidates = [(i, j) for i in range(k, n) for j in range(k, n)
                          if not _zero_possible(a[i][j])]
            if not candidates:
                if local:
                    raise PrecisionError(f"no certified nonzero pivot at step {k}; "
                                         "singularity is unresolved at this precision")
                raise ZeroDivisionError("matrix is singular")
            i, j = min(candidates, key=lambda ij: a[ij[0]][ij[1]].v) if local else candidates[0]
            if local and any(a[r][c].v < a[i][j].v for r in range(k, n) for c in range(k, n)):
                raise PrecisionError(f"cannot certify a minimum-valuation pivot at step {k}; "
                                     "a zero-containing ball has a lower valuation bound")
            a[k], a[i] = a[i], a[k]
            b[k], b[i] = b[i], b[k]
            row_order[k], row_order[i] = row_order[i], row_order[k]
            for row in a:
                row[k], row[j] = row[j], row[k]
            column_order[k], column_order[j] = column_order[j], column_order[k]
            pivot = a[k][k]
            if local:
                pivots.append(pivot.valuation())
            for r in range(k+1, n):
                factor = a[r][k]/pivot
                for c in range(k+1, n):
                    a[r][c] = a[r][c]-_multiply(factor, a[k][c])
                # This zero is a structural consequence of the row operation.
                a[r][k] = Fraction(0)
                for c in range(width):
                    b[r][c] = b[r][c]-_multiply(factor, b[k][c])
        x = [[Fraction(0)]*width for _ in range(n)]
        for i in range(n-1, -1, -1):
            for c in range(width):
                numerator = b[i][c]-_dot(a[i][i+1:], [x[j][c] for j in range(i+1, n)])
                x[i][c] = Fraction(0) if _exact_zero(numerator) else numerator/a[i][i]
        ordered = [None]*n
        for i, original in enumerate(column_order):
            ordered[original] = x[i]
        if local:
            # Exact zero solutions may be enclosed at any precision. Do not
            # cap them below the known input/output digits at a field default.
            zero_precision = max([self.base.precision]
                                 + [c.precision for row in self.rows for c in row]
                                 + [c.precision for row in ordered for c in row if hasattr(c, "precision")])
            ordered = [[self.base(0, absprec=zero_precision) if _exact_zero(c) else self.base(c)
                        for c in row] for row in ordered]
        solution = tuple(row[0] for row in ordered) if vector_rhs else Matrix(self.base, ordered)
        precisions = (tuple(getattr(c, "precision", None) for c in solution) if vector_rhs
                      else solution.precisions)
        e = (1 if isinstance(self.base, PadicField) else self.base.ramification_index) if local else None
        return LinearSolveResult(solution, tuple(row_order), tuple(column_order), tuple(pivots), precisions, e)
