Precision model and guarantees
==============================

A numerical value denotes a ball

.. math::

   x = a + p^N\mathbb Z_p.

The stored ``precision`` is N, ``v`` is v_p(a) for a known nonzero leading
digit, and ``relative_precision`` is N-v. For a zero ball v=N is only a
lower bound. This follows the usual absolute/relative precision distinction
described in the `Sage p-adics tutorial <https://doc.sagemath.org/html/en/reference/padics/sage/rings/padics/tutorial.html>`_.

Arithmetic enclosures
---------------------

For balls with precision N_x, N_y and valuation lower bounds v_x, v_y:

.. math::

   N_{x+y}=\min(N_x,N_y),\qquad
   N_{xy}=\min(N_x+v_y,N_y+v_x).

For a ball excluding zero:

.. math::

   v_{1/x}=-v_x,\qquad N_{1/x}=N_x-2v_x.

These follow by expanding perturbations, including their product. The
inverse bound follows from ``1/(a+h)-1/a=-h/(a*(a+h))``. Exact rational
scalar multiplication shifts N by the scalar valuation; division shifts it
in the opposite direction. No float enters these computations.

The enclosures are conservative. The implementation does not track correlated
errors: for example, ``x-x`` produces O(p^N), not an exact symbolic zero.
Repeated squaring can lose information that a specialized differential
precision algorithm might retain. A new exact source can be evaluated at
higher precision, but an existing finite ball cannot be extended by choosing
its missing digits to be zero.

Analytic bounds
---------------

Exponential, sine and cosine use

.. math::

   v_p(x^k/k!)\ge k\left(v_p(x)-\frac1{p-1}\right).

Once this lower bound reaches the target N, every later omitted term does
too. In their supported domains the output is a safe enclosure at input
absolute precision. The logarithm on ``1+p Zp`` uses

.. math::

   v_p(x^k/k)\ge k v_p(x)-\lfloor\log_p k\rfloor.

The right-hand side is nondecreasing for integral v_p(x)>=1. A small isolated
term is therefore never used as evidence that the entire tail is small.
For general nonzero x, removing its p-power shifts log's input precision
from N to N-v_p(x).

For square roots, the derivative gives output precision N-v_p(root) for
odd p, and N-v_p(root)-1 at p=2. The implementation checks the square-unit
conditions before lifting. See :doc:`functions` for theta and zeta bounds.

Limits of the guarantee
-----------------------

The precision formulas are mathematical enclosures implemented and tested
with exact integer/rational arithmetic; the code is not formally verified.
Complex images use mpmath numerical precision and are not certified complex
intervals. Quadratic extensions use coefficient balls, not a full
uniformizer-based extension-field precision model. Exact rational algebra
and numerical balls are distinct layers.
