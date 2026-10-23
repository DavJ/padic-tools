"""Exact and certified local-field linear systems; no optional dependencies."""
from fractions import Fraction

from padic_tools import Gaussian, GaussianLocalPlace, Matrix, Qp, exact


def main():
    source = Matrix(exact, [[1, 2], [3, 7]])
    assert source.solve([1, 0]) == (Fraction(7), Fraction(-3))
    print("Exact rational solution:", source.solve([1, 0]))

    conditioned = Matrix(Qp(3, 12), [[9, 0], [0, 243]], absprec=21)
    report = conditioned.solve_with_info([1, 1])
    assert all(x.contains(y) for x, y in zip(report.solution, (Fraction(1, 9), Fraction(1, 243))))
    print("Q3 pivot valuations:", report.pivot_valuations)
    print("Q3 condition exponent:", report.condition_exponent)
    print("Q3 solution precisions:", report.precisions)
    print("Q3 inverse precisions:", conditioned.inverse().precisions)

    gaussian = Matrix(Gaussian, [[Gaussian(1, 1), 2], [3, Gaussian(2, -1)]])
    place = GaussianLocalPlace(Qp(2, 12), root_sign=-1)
    mapper = lambda c: place(c, absprec=23)
    ramified = gaussian.map_coefficients(place.field, mapper)
    assert ramified.inverse().is_congruent(gaussian.inverse().map_coefficients(place.field, mapper))
    assert ramified.determinant().is_congruent(mapper(gaussian.determinant()))
    print("Q2(i) inverse agrees with the mapped exact source.")


if __name__ == "__main__":
    main()
