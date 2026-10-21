# Validation of the 0.1 implementation

Validated on 22 October 2026 with Python 3.12.14.

- **121 pytest cases passed**, including all seven independent PARI/GP
  comparisons. GP 2.15.2 was run from a locally extracted Debian package;
  the package SHA-256 was checked against Debian's package index. The native
  evaluator agreed with PARI's principal branch at matched negative and
  positive integer arguments for p=2,3,5,7.
- Exact rational-oracle checks cover arithmetic and perturbations of ball
  centres, inverse precision, cancellation, valuations and negative absolute
  precision. Analytic checks include exp/log identities, independently
  summed rational exponential series, dyadic boundaries and square roots.
- Forty zeta branch/negative-integer combinations agree with exact
  Bernoulli/Euler-factor values, including nonprincipal even branches.
- Bridge tests check Gaussian and quaternion sums/products/inverses across
  coefficient base changes, chosen split places, zero divisors and the
  explicit even-Clifford isomorphism. Clifford associativity is tested on
  all triples of basis blades in a mixed-signature example.
- **84 documentation examples passed** using Sphinx doctest. HTML builds
  with warnings treated as errors; local links/assets were checked, and all
  12 displayed equations were parsed to MathML by the bundled KaTeX renderer.
- Wheel and source archive built successfully, and **twine check passed**.
  The wheel installed and its core/bridge/zeta smoke checks worked in a fresh
  virtual environment with **no optional dependencies**.

Run the commands in README and RELEASING to reproduce. The optional PARI
tests skip when `gp` is not installed; CI installs it explicitly. The
Python-version matrix is configured in GitHub Actions; this local record
does not claim those remote jobs have already passed.

These checks support the documented 0.1 scope; they do not establish formal
verification, a general complex/p-adic analytic comparison theorem, a full
extension-field engine, or production performance guarantees.

## Portable documentation repair

The documentation now includes a reproducible single-file HTML edition with
embedded styles and pre-rendered MathML. Its builder validates all **222
internal anchors**, removes controls requiring companion pages, and checks the
English/Czech guide targets. CI builds and uploads this edition separately.

Validated in Chromium 150 with network access and JavaScript disabled at
1440x1000 and 390x844: English/Czech navigation, return to the English start,
and chapter navigation worked; all 12 formulas were native MathML. The normal
multi-page edition's English/Czech links also worked from its complete local
directory. Both Sphinx HTML builds completed without warnings and all 84
documentation examples passed. The source distribution includes the custom
language-navigation template and portable manual builder.
