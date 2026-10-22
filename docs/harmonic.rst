Harmonic analysis on Qp
=======================

The harmonic layer integrates functions of a p-adic variable with values
in a chosen complex place. Its native source coefficients are rational
numbers and exact rational combinations of p-power roots of unity.
``PadicTestFunction`` specifies the complete data of a compactly supported,
locally constant function. All Fourier sums, integrals and identities below
are exact; optional complex evaluation is a separate numerical operation.

Characters and known digits
---------------------------

Our additive character has positive phase and kernel Zp:

.. math::

   \chi_p(x)=\exp(2\pi i\{x\}_p),\qquad
   \chi_p(x+y)=\chi_p(x)\chi_p(y),\qquad
   \chi_p(x)=1\ \Longleftrightarrow\ x\in\mathbb Z_p.

``additive_phase`` returns the unique rational q in [0,1), with p-power
denominator, such that x-q belongs to Zp. This uses the negative p-adic
digits. A rational's real fractional part generally gives a different
answer. A ``Padic`` ball must have absolute precision at least zero to
determine this phase; less information raises ``PrecisionError``.

.. doctest::

   >>> K = Qp(3, 18)
   >>> additive_phase(Fraction(1, 6), K)
   Fraction(2, 3)
   >>> eta = additive_character(Fraction(1, 9), K)
   >>> eta**9 == 1
   True
   >>> eta*eta.conjugate() == 1
   True
   >>> additive_character(K(7)) == 1
   True
   >>> additive_character(K.O(-1))
   Traceback (most recent call last):
       ...
   padic_tools.core.PrecisionError: additive character requires absolute precision >=0

The output of ``additive_character`` is an exact ``CyclotomicValue``.
It preserves addition of arguments as multiplication of characters.
This map is not a field embedding into C.

Exact cyclotomic coefficients
-----------------------------

``CyclotomicValue(p, coefficients, level=n)`` represents a polynomial in
the specified complex root zeta=exp(2*pi*i/p**n), with rational coefficients
in ascending order. The native algorithm reduces by

.. math::

   \Phi_{p^n}(X)=\sum_{j=0}^{p-1}X^{j p^{n-1}},\qquad n\geq1.

It then descends to the smallest p-power root field containing the value.
Consequently equality, zero sums, conjugation and rational scaling are
exact, including when two formulas start at different root levels. The
chosen complex embedding fixes the phase; no p-adic root place is inferred.
Native operations cover the coefficient ring, rather than a general
number-field API with algebraic inverses.

.. doctest::

   >>> i = CyclotomicValue.root(2, 2)
   >>> i*i == -1
   True
   >>> sum(CyclotomicValue.root(3, 2, j) for j in range(9)) == 0
   True
   >>> i.to_sympy()
   I
   >>> abs(i.complex(dps=60).imag-1) < 10**-55
   True

``.to_sympy()`` requires the symbolic extra. ``.complex(dps=...)`` uses
the complex extra and an isolated mpmath context. Its dps is a working
precision, not a certified complex error enclosure. Exact source data is
retained; floats, Python complex values and p-adic balls are rejected as
function coefficients.

Real Haar integration
----------------------

We normalize real additive Haar measure by vol(Zp)=1. Translation and
dilation give the exact volume of any ball:

.. math::

   \operatorname{vol}(c+p^b\mathbb Z_p)=p^{-b}.

This real measure underlies the Fourier theory in this chapter. Its weights
are distinct from the Qp-valued bounded measures of :doc:`measures` and the
Mellin/Amice integrals in :doc:`transforms`.

.. doctest::

   >>> haar_volume(K, 2), haar_volume(K, -2)
   (Fraction(1, 9), Fraction(9, 1))
   >>> unit_ball = PadicTestFunction.ball(K)
   >>> unit_ball.integral().as_fraction()
   Fraction(1, 1)
   >>> small_ball = PadicTestFunction.ball(K, center=Fraction(1, 3), exponent=1)
   >>> small_ball(Fraction(1, 3)) == 1 and small_ball(0) == 0
   True
   >>> small_ball.integral().as_fraction()
   Fraction(1, 3)

A complete finite model on Qp
-----------------------------

Let a<=b be integers. A function supported inside p**a*Zp and constant on
cosets of p**b*Zp is completely described by N=p**(b-a) values:

.. math::

   f(x)=f_j\quad(x\in p^a j+p^b\mathbb Z_p),\qquad
   j=0,\ldots,N-1,\qquad
   \int f(x)\,dx=p^{-b}\sum_{j=0}^{N-1}f_j.

These exponents specify a support container and a constancy lattice, not a
numerical cutoff or an optimal support description. Values outside the
container are zero. Point evaluation accepts exact rational data or a ball
whose known digits determine the relevant coset. A ball known to lie
outside the container returns zero immediately. Construction needs every
coset value; it never fills an unspecified tail with zeros.

.. doctest::

   >>> f = PadicTestFunction(K, [1, 2, 4], support_exponent=0, constancy_exponent=1)
   >>> f(1) == 2 and f(4) == 2 and f(Fraction(1, 3)) == 0
   True
   >>> f.integral().as_fraction(), f.norm_squared().as_fraction()
   (Fraction(7, 3), Fraction(7, 1))

``field.precision`` does not round these exact source values. Operations
on two test functions, including pointwise products, inner products and
convolution, require the same prime and the same two lattice exponents.
``translate(t)`` represents f(x-t) for t inside the support container;
``reflect()`` represents f(-x). ``on_lattice`` exactly enlarges the container
and/or refines the cosets; coarsening is rejected. Subtraction and negation
are supported on a common lattice.

Fourier transform and its normalization
---------------------------------------

Our Fourier kernel is negative, and the inverse kernel is positive:

.. math::

   \widehat f(\xi)=\int_{\mathbb Q_p}f(x)\chi_p(-x\xi)\,dx.

Periodicity under p**b*Zp forces the transform to vanish outside
p**(-b)*Zp: translating by an element of that subgroup multiplies the
integral by a character, so a nontrivial factor forces zero. Support
inside p**a*Zp makes the transform constant modulo p**(-a)*Zp.
On the resulting dual grid the integral is the finite sum

.. math::

   \widehat f(p^{-b}k)=p^{-b}\sum_{j=0}^{N-1}f_j
      \exp\!\left(-2\pi i\frac{jk}{N}\right),\qquad
   (a,b)\longmapsto(-b,-a).

Each source coset has Haar volume p**(-b), and its character is constant
at these dual frequencies; integrating over that coset proves the formula.
There is no numerical integration error for this declared function class.

.. doctest::

   >>> transformed = f.fourier()
   >>> transformed.support_exponent, transformed.constancy_exponent
   (-1, 0)
   >>> transformed(0) == f.integral()
   True
   >>> transformed.fourier(inverse=True) == f
   True
   >>> transformed.fourier() == f.reflect()
   True
   >>> unit_ball.fourier() == unit_ball
   True

For raw cyclic arrays, ``finite_fourier(values, K)`` instead uses forward
normalization 1/N and an unnormalized inverse. Its kernel has the same
sign. ``PadicTestFunction.fourier`` supplies the Haar scaling and dual-grid
metadata needed for the full Qp integral. These two normalization contracts
are intentionally explicit.

Parseval and convolution
-------------------------

The finite geometric sum is N when its frequency index is divisible by N,
and zero otherwise. Applying it to the two finite sums proves inversion
and the Haar Parseval identity; multiplying and regrouping proves the
convolution theorem:

.. math::

   \int f(x)\overline{g(x)}\,dx
     =\int\widehat f(\xi)\overline{\widehat g(\xi)}\,d\xi,
   \qquad
   \widehat{f*g}=\widehat f\,\widehat g,\qquad
   (f*g)_k=p^{-b}\sum_{j=0}^{N-1}f_jg_{k-j\bmod N}.

The dual-coset volume is p**a. Together with N=p**(b-a), this cancels
the squared p**(-b) Fourier weight in Parseval. Thus the formula holds
also for negative exponents and support larger than Zp.

.. doctest::

   >>> transformed.norm_squared() == f.norm_squared()
   True
   >>> g = PadicTestFunction(K, [3, 0, 1])
   >>> f.inner_product(g) == f.fourier().inner_product(g.fourier())
   True
   >>> f.convolve(g).fourier() == f.fourier()*g.fourier()
   True
   >>> spectrum = finite_fourier([1, 2, 4], K)
   >>> finite_fourier(spectrum, K, inverse=True) == (1, 2, 4)
   True

Work limits, validation and next steps
---------------------------------------

Test-function construction defaults to at most 4096 cosets. The native
direct Fourier and convolution methods default to at most 256 cosets;
``max_size`` can be raised explicitly. The exact coefficient backend has
a fixed maximum root order of 4096. Resource limits raise ``PrecisionError``
and do not silently discard modes. The direct DFT uses quadratic work for
rational input; dense cyclotomic source values add coefficient work.

Tests check exact inversion, both Parseval identities, the translation and
convolution theorems, and the analytic Fourier transform of ball indicators.
Independent SymPy polynomial quotients check cyclotomic reduction, while
PARI/GP checks the complex embedding of finite spectra. Character, kernel
and Haar conventions follow Tao's p-adic discussion; Kedlaya gives the
general duality and inversion context. Both sources are in :doc:`references`.

This foundation supports finite prime-sector spectra and wave packets.
Kozyrev bases, whole-Qp Vladimirov images and conservative ball diffusion
are described in :doc:`wavelets`. Multidimensional grids and general
complex sample coefficients remain subsequent extensions. Applying the
API to a function outside its declared class requires a separate error
estimate for the chosen approximation.
