# padic-tools

A Python library for finite precision **p-adic analysis**, exact **symbolic
algebra**, and explicit bridges to **complex numbers, quaternions,
biquaternions and Clifford algebras**.

This is a usable **0.1 research library**, not a claim to implement every part
of p-adic mathematics. Domains, precision bounds and remaining scope are
documented. The numerical core uses only Python's standard library.

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
| Analytic | exp, Iwasawa log with explicit log(p), log1p, sin, cos, principal-unit powers, square roots including p=2 |
| Special functions | Multiplicative theta3, theta4 constants, exact Bernoulli/zeta values, numerical even KL zeta branches on Zp away from s=1 |
| Symbolic | Rational polynomials, Euclidean division/gcd, formal series, composition/calculus, optional SymPy function nodes/evaluation |
| Algebraic bridges | Q and Q[i] to independent completions, chosen split Gaussian places, rational/p-adic quaternions and biquaternions |
| Geometric algebra | Diagonal Clifford algebras, grades, geometric/wedge products, involutions, explicit Hamilton-even-Clifford isomorphism |

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

