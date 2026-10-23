# padiclab

A Python library for finite precision **p-adic analysis**, exact **symbolic
algebra**, and explicit bridges to **complex numbers, quaternions,
biquaternions and Clifford algebras**.

This is a usable **0.1 research library**, not a claim to implement every part
of p-adic mathematics. Domains, precision bounds and remaining scope are
documented. The numerical core uses only Python's standard library.

Repository: [DavJ/padiclab](https://github.com/DavJ/padiclab).
The current Python distribution is `padic-tools`, imported as `padic_tools`.

```bash
python -m pip install .
# Optional symbolic expressions, complex evaluation and HTML documentation:
python -m pip install '.[symbolic,complex,docs]'
```

```python
from fractions import Fraction
from padic_tools import Qp, exp, log, theta3, zeta_p

K = Qp(5, precision=20)
x = K(Fraction(1, 3))
assert (x * 3).is_congruent(1)
assert log(exp(K(5))).is_congruent(5)

t = theta3(K(5))
z = zeta_p(-3, K, branch=0)  # Euler-corrected zeta(-3), on this branch
assert z.is_congruent(Fraction(-31, 30))
```

| Layer | Implemented in 0.1 |
|---|---|
| Numerical | Qp balls, valuations, digits, inverses, precision propagation, simple-root Hensel lifting |
| Local fields | Certified unramified and Eisenstein quadratic extensions, Q3(i) and Q2(i), ideal-ball precision, exact residue fields and norm/trace; lifted Frobenius for unramified fields |
| Analytic | exp, Iwasawa log with explicit log(p), log1p, sin, cos, principal-unit powers, square roots including p=2 |
| Special functions | Multiplicative theta3, theta4 constants, exact Bernoulli/zeta values, numerical even KL zeta branches on Zp away from s=1 |
| Symbolic | Rational polynomials, Euclidean division/gcd, formal series, composition/calculus, optional SymPy function nodes/evaluation |
| Transforms | Exact atomic and bounded rational Amice measures, moments/polynomial integration, unit restriction, certified branched Mellin, local Laplace, Amice evaluation and additive convolution |
| Taylor | Exact polynomial shifts and rational jets; local p-adic exp/log/log1p/sin/cos and polynomial expansions with whole-tail enclosures |
| Harmonic analysis | Exact additive characters, real Haar integration, finite Fourier transforms, compactly supported test functions, inversion, Parseval and convolution |
| Wavelets and operators | Exact finite Kozyrev bases and reconstruction, whole-Qp integer-order Vladimirov images with analytic tails, conservative ball generators and numerical diffusion |
| Algebraic bridges | Q and Q[i] to independent completions; chosen split, unramified and ramified local-field images; rational quaternions and biquaternions |
| Geometric algebra | Diagonal Clifford algebras, grades, geometric/wedge products, involutions, explicit Hamilton-even-Clifford isomorphism |

## Unramified quadratic fields

```python
from padic_tools import UnramifiedQuadraticField

L = UnramifiedQuadraticField(Qp(3, 12), [1, 0, 1], name="i")
x = L(1, 2)
assert (L.generator**2).is_congruent(-1)
assert (x*x.inverse()).is_congruent(1)
assert x.norm().contains(5) and x.trace().contains(2)
assert L.residue_field().order == 9
assert x.frobenius().residue() == x.residue()**3
```

The constructor checks monicity, p-integral exact coefficients and
irreducibility modulo p. `[1, 0, 1]` defines Q3(i); it is rejected over Q5
(split) and Q2 (ramified). `[1, 1, 1]` instead defines the unramified
quadratic extension of Q2. Each element is a whole ideal ball with common
coordinate precision, including negative valuations. `norm()` is algebraic;
`absolute_value()` is the p-adic magnitude. See [the extension chapter](docs/extensions.rst)
for error bounds, exact residues and scope. Higher degrees, extension
analytic functions are subsequent steps; dedicated field adapters are included.

## Ramified Q2(i) and Eisenstein fields

```python
from padic_tools import Q2i, EisensteinQuadraticField

L2i = Q2i(precision=12)
pi = L2i.uniformizer()  # pi = 1+i, with pi^2-2*pi+2 = 0
i = pi-1
assert (i*i).is_congruent(-1)
assert pi.valuation() == 1 and L2i(2).valuation() == 2
assert pi.base_valuation() == Fraction(1, 2)
assert pi.norm().contains(2) and pi.trace().contains(2)
assert L2i.residue_field().order == 2

ball = L2i(1, 1, absprec=5)  # Precision is in powers of pi
assert tuple(c.precision for c in ball.coefficients) == (3, 2)
R = EisensteinQuadraticField(Qp(3, 12), [-3, 0, 1])
assert (R.generator**2).is_congruent(3)
```

Coordinates use the basis `(1, pi)`. Odd and negative uniformizer precisions
give staggered base-coordinate precisions. Norm and trace return Qp balls
in base p units, with bounds that account for the different, including wild
dyadic ramification. `absolute_value_squared()` gives the exact rational
square of the base-compatible magnitude. See [the ramified chapter](docs/ramified.rst)
for the derivation and explicit Gaussian source conversion.

## Chosen Gaussian and quaternion local places

```python
from padic_tools import Gaussian, GaussianBridge, GaussianLocalPlace, Biquaternions
from padic_tools import QuaternionBridge, CliffordAlgebra

place = GaussianLocalPlace(Qp(2, 12), root_sign=1)  # i -> pi-1 in Q2(i)
source = GaussianBridge(Fraction(1, 3), 2)
image = source.local(place, absprec=19)  # Absolute powers of pi
assert image.precision == 19 and image.contains(Fraction(-5, 3), 2)
assert image.norm().contains(source.value.norm())

q = Biquaternions(Gaussian(1, 2), Gaussian(3, -1), 0, 1)
local_q = QuaternionBridge(q).local(place)
C = CliffordAlgebra((1, 1, 1), place.field)
assert q.to_clifford().map_coefficients(C, place).is_congruent(local_q.to_clifford())
```

The same API selects Qp at split primes and an unramified quadratic field
at primes congruent to 3 modulo 4. `root_sign` explicitly chooses the
generator image. See [chosen local places](docs/local_bridges.rst) for
guard digits, conjugations, algebraic norms and the Clifford commuting square.
Extension analytic functions remain a later step.

## Transforms and Taylor expansions

```python
from padic_tools import AtomicMeasure, RationalAmiceMeasure, taylor

mu = AtomicMeasure(K, [(1, 1), (2, 2)])
assert mu.mellin(6, branch=2).is_congruent(mu.moment(6))
assert mu.amice(exp(K(-5))-1).is_congruent(mu.laplace(K(5)))

jet = taylor(exp, K(5), order=30)
assert jet(K(10)).is_congruent(exp(K(10)), precision=20)

rational = RationalAmiceMeasure(K, [1], [2, 1])  # Complete A(T)=1/(2+T)
assert rational.moment(3) == Fraction(1, 8)
assert rational.laplace(K(5)).is_congruent(1/(1+exp(K(-5))))

units = rational.restrict_units()  # Exact rational transform of the restriction
assert units.moment(3) == Fraction(-31, 2)
assert units.mellin(3, branch=3).contains(units.moment(3))
assert units.mellin(Fraction(1, 2), branch=3).precision == 20
```

The measure and kernel convention are explicit. Mellin support must lie in
Zp units; Laplace checks the local p-adic exponential domain. Atomic measures
are finite exact data. Rational Amice measures additionally describe complete
bounded moment sequences with a checked coefficient bound. Their Mellin
evaluator proves an infinite Mahler-tail bound, including the dyadic case;
`precision` requests an absolute output order and `max_terms` limits work.
General bounded-measure oracles and transform inverses are future extensions.
Formal series order and p-adic numerical precision are
tracked separately; a numerical Taylor expansion includes its omitted tail.

## Exact harmonic analysis

```python
from padic_tools import PadicTestFunction, additive_phase

K3 = Qp(3)
assert additive_phase(Fraction(1, 6), K3) == Fraction(2, 3)
f = PadicTestFunction(K3, [1, 2, 4])  # Supported on Z3, constant modulo 3
assert f.integral().as_fraction() == Fraction(7, 3)
spectrum = f.fourier()
assert spectrum.fourier(inverse=True) == f
assert spectrum.norm_squared() == f.norm_squared()
assert PadicTestFunction.ball(K3).fourier() == PadicTestFunction.ball(K3)
```

The coefficients are exact rational combinations of complex p-power roots of
unity; `.complex(dps=...)` provides an optional numerical view. Haar weights
are real, with volume(Zp)=1. The transform retains support/constancy metadata
for the dual Qp lattice. Raw `finite_fourier` arrays have an explicit 1/N
forward normalization. Native work limits are documented in the harmonic
chapter; accelerated backends and general complex samples remain future work.

## Wavelets and conservative diffusion

```python
from padic_tools import KozyrevBasis, kozyrev_wavelet, vladimirov, BallVladimirov

basis = KozyrevBasis(K3, constancy_exponent=2)
signal = PadicTestFunction(K3, range(9), constancy_exponent=2)
assert basis.reconstruct(basis.decompose(signal)) == signal
wave = kozyrev_wavelet(K3, scale=1)
assert vladimirov(wave, order=2).as_test_function() == wave * 81

operator = BallVladimirov(K3, constancy_exponent=2)
initial = (PadicTestFunction.ball(K3, exponent=2)*9).on_lattice(
    support_exponent=0, constancy_exponent=2)
state = operator.evolve(initial, Fraction(1, 10), dps=60)  # complex extra
assert abs(state.integral()-1) < 10**-50
```

Wavelets retain their exact squared norms. `vladimirov` includes the entire
exterior tail for nonzero-mass input on Qp. `BallVladimirov` instead uses the
regional conservative boundary on its declared ball. Its exact spectrum and
kernel matrix determine mass-preserving heat evolution. Numerical outputs
use isolated mpmath working precision. Run `python examples/padic_diffusion.py`
for a reproducible positivity, mass, energy and grid/boundary comparison.

## The bridge is included

```python
from padic_tools import Gaussian, GaussianBridge, Hamilton, QuaternionBridge

g = GaussianBridge(Gaussian(Fraction(1, 3), 2))
complex_image = g.complex(dps=60)
local_algebra_image = g.padic_algebra(Qp(3, 20))
chosen_split_image = g.split_padic(Qp(5, 20), root_sign=1)

q = Hamilton(1, 2, 3, 4)
q5 = QuaternionBridge(q).padic(Qp(5, 20))
assert (q*q).to_clifford() == q.to_clifford()*q.to_clifford()
```

The common source is an exact algebra over **Q**; its complex and p-adic
images are calculated separately. A finite p-adic residue is not a complex
approximation. At split primes the quadratic quotient has zero divisors.
The complex Riemann zeta and a KL branch are distinct functions with a
specified special-value interpolation relation. See the bridge chapter for
proofs, conventions and examples.

## Documentation and validation

```bash
python -m pip install '.[dev]'
python -m pytest
python -m sphinx -W --keep-going -b html docs docs/_build/html
python -m sphinx -W --keep-going -b doctest docs docs/_build/doctest
```

Open **docs/_build/html/index.html**. The guide includes a Czech introduction,
an API reference, mathematical derivations, executable examples, a feature
roadmap and links to primary sources. CI uploads the generated HTML manual as
the `html-documentation` artifact; it does not publish a website.

For sharing or opening a single file, build the portable manual instead:

```bash
python scripts/build_manual.py
```

Open **docs/_build/padic-tools-manual.html**. This file includes every chapter,
its styles and pre-rendered mathematical formulas. Chapter and language links
stay inside the file, so it works offline without companion files or
JavaScript. Building this portable edition additionally requires Node.js;
reading it does not. CI uploads it as `portable-html-manual`.

The regular multi-page HTML edition must be served as a complete directory
or opened from the full downloaded `html-documentation` archive. Sharing
only its `index.html` loses the other pages and assets. Both editions provide
a separate **English · Čeština** language navigation bar; the Czech guide
is kept out of the English chapter menu.

Start with [the Czech guide](docs/cesky.rst),
[the numerical tutorial](docs/tutorial.rst),
[the special-function definitions](docs/functions.rst), and
[the complex/quaternion/Clifford bridge](docs/bridges.rst).

The 2022 `main.py` helpers remain as an installed-source compatibility demo.
`padic()` now returns a ball; new code should import `padic_tools`.
The unfinished float-based arithmetic has been replaced with exact rational
conversion. The project does not currently grant an open-source license;
licensing remains the repository owner's decision.

