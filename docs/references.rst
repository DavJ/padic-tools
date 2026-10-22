Sources and conventions
========================

Primary documentation and research sources used to define and check the
implementation are linked below. The bridge identification and tail-bound
derivations are written explicitly in this guide. No external source is
presented as establishing a new complex/p-adic zero correspondence.

* `Sage: Introduction to p-adics <https://doc.sagemath.org/html/en/reference/padics/sage/rings/padics/tutorial.html>`_
  — absolute versus relative precision, valuation and numerical representations.
* `Sage: p-adic generic elements <https://doc.sagemath.org/html/en/reference/padics/sage/rings/padics/padic_generic_element.html>`_
  — analytic domains, logarithm branches and precision behavior.
* `Sage: Unramified Extension Generic
  <https://doc.sagemath.org/html/en/reference/padics/sage/rings/padics/unramified_extension_generic.html>`_
  — residue fields, generators, Galois structure and the unchanged uniformizer
  in unramified extensions. The quadratic valuation and error derivations
  for our ideal balls are given in :doc:`extensions`.
* `Magma: Overview of the p-adics
  <https://docs.magma-maths.org/LocalFields/PAdicRingsAndTheirExtensions/overview.html>`_
  — absolute precision in uniformizer units and its scaling by the
  ramification index when moving between local fields.
* `NIST DLMF 20.2 <https://dlmf.nist.gov/20.2>`_
  — Jacobi theta definitions; our multiplicative argument convention is
  distinguished in :doc:`functions`.
* `NIST DLMF 25.6 <https://dlmf.nist.gov/25.6>`_
  — Bernoulli formulas for classical zeta special values.
* `Rodrigues Jacinto and Williams: An Introduction to p-adic L-functions
  <https://warwick.ac.uk/fac/sci/maths/people/staff/cwilliams/lecturenotes/lecture_notes_part_i.pdf>`_
  — the KL construction, special-value interpolation and character variables;
  see Theorem 3.12 and section 4.3.
  Sections 2.2--2.4 also describe the measure/series correspondence,
  binomial moments, Dirac measures and the root-of-unity restriction formula
  used in :doc:`transforms` and :doc:`mellin`.
* `Roblot: Computing p-adic L-functions of totally real number fields
  <https://arxiv.org/html/1110.0246v1>`_
  — Lemma 5.8, Proposition 5.10 and Corollary 5.13 give Mahler integration
  and kernel tail bounds. The adaptation to every torsion branch, and the
  implementation's precision contract, are explained in :doc:`mellin`.
* `Tao: Tate's proof of the functional equation
  <https://terrytao.wordpress.com/2008/07/27/tates-proof-of-the-functional-equation/>`_
  — the p-adic section defines normalized real Haar measure, the positive
  additive character, the negative Fourier kernel and the self-dual unit-ball
  indicator. These conventions agree with :doc:`harmonic`.
* `Kedlaya: Adelic Fourier analysis after Tate
  <https://kskedlaya.org/cft/sec_adelic-fourier.html>`_
  — section 6.6 describes additive duality and local Fourier inversion,
  including the dependence on character and Haar normalization. Our finite
  lattice derivation is given in :doc:`harmonic`.
* `Sage: Algebraic numbers and number fields
  <https://doc.sagemath.org/html/en/reference/number_fields/sage/rings/number_field/number_field.html>`_
  — cyclotomic fields and their chosen complex embeddings. The native
  harmonic coefficient ring implements rational combinations of p-power
  roots; it does not provide general number-field arithmetic.
* `Kozyrev: Wavelet analysis as a p-adic spectral analysis
  <https://arxiv.org/html/math-ph/0012019v3>`_
  — the wavelet eigenbasis and Vladimirov kernel/multiplier convention.
  Our opposite Fourier sign, reversed scale and amplitude-one normalization
  are stated in :doc:`wavelets`.
* `Bikulov and Zubarev: On one real basis for L2(Qp)
  <https://arxiv.org/abs/1504.03624>`_
  — compact-ball Fourier analysis and the regional integral, distinct from
  the whole-space derivative of a zero extension. Our sign, kernel constant
  and escape-rate derivation appear in :doc:`wavelets`.
* `Mathlib: The Amice transform
  <https://leanprover-community.github.io/mathlib4_docs/Mathlib/NumberTheory/Padics/Measure/AmiceTransform.html>`_
  — the formal definition and coefficient identity for the Amice transform;
  consulting this reference is not a claim that our Python algorithms have
  been formally verified in Lean.
* `Lai, Lupu and Sprang: On the irrationality of certain p-adic zeta values
  <https://link.springer.com/article/10.1007/s40687-025-00559-x>`_
  — sections 2.1–2.3 recall Volkenborn moments and Teichmüller/Hurwitz
  definitions. Their integer-argument zeta convention is not silently
  identified with a fixed ``branch=0`` for all s.
* `PARI/GP: Transcendental functions
  <https://pari.math.u-bordeaux.fr/dochtml/html-stable/Transcendental_functions.html>`_
  — an independent p-adic zeta implementation; input and branch conventions
  must be matched before numerical comparison.
* `Sage: Quaternion algebras
  <https://doc.sagemath.org/html/en/reference/algebras/sage/algebras/quatalg/quaternion_algebra.html>`_
  — quaternion algebra construction over a base field.
* `Sage: Clifford algebras
  <https://doc.sagemath.org/html/en/reference/algebras/sage/algebras/clifford_algebra.html>`_
  — general Clifford relations, dimension and universal property.
* `Sage: Number-field embeddings
  <https://doc.sagemath.org/html/en/reference/number_fields/sage/rings/number_field/number_field_morphisms.html>`_
  — explicit embeddings into ambient coefficient domains.
* `Deterministic Miller–Rabin base records <https://miller-rabin.appspot.com/>`_
  — the seven test bases used for integers below 2^64; larger inputs are
  rejected by the native prime validator.
* `SymPy: Writing custom functions <https://docs.sympy.org/latest/guides/custom-functions.html>`_
  and `mpmath contexts <https://mpmath.org/doc/current/contexts.html>`_
  — optional symbolic nodes and isolated archimedean precision contexts.
* `PyPI: Publishing to PyPI <https://docs.pypi.org/publishing/publishing-to-pypi/>`_
  and `Trusted Publishing <https://docs.pypi.org/trusted-publishers/using-a-publisher/>`_
  — packaging and the optional manual release workflow.
