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
