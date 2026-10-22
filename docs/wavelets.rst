Wavelets, Vladimirov operators and diffusion
============================================

This layer builds on the exact Haar/Fourier conventions in :doc:`harmonic`.
Wavelets, decomposition coefficients and positive integer-order operators
use rational/cyclotomic arithmetic. Heat evolution uses a separate complex
numerical view with an explicit working precision.

Kozyrev wavelets and finite bases
---------------------------------

For a scale c, centre z and frequency j=1,...,p-1, the API uses

.. math::

   w_{c,z,j}(x)=\chi_p\!\left(jp^{-c-1}(x-z)\right)
       1_{z+p^c\mathbb Z_p}(x),\qquad
   \int w_{c,z,j}=0,\qquad \|w_{c,z,j}\|_2^2=p^{-c}.

These are amplitude-one versions of Kozyrev's wavelets. Our c equals
minus his gamma; centring the character changes a constant unit phase.
We retain exact squared norms instead of adjoining square roots for unit
normalization. The coefficient ring therefore remains the one declared
in :doc:`harmonic`. A ball centre must be known through order c+1 to
determine both its support and its phase.

.. doctest::

   >>> K = Qp(3, 18)
   >>> wave = kozyrev_wavelet(K, scale=1, frequency=1)
   >>> wave.integral().as_fraction(), wave.norm_squared().as_fraction()
   (Fraction(0, 1), Fraction(1, 3))
   >>> wave(0) == 1 and wave(1) == 0
   True

On a lattice with support container p**a*Zp and constancy p**b*Zp,
``KozyrevBasis`` contains the container indicator plus all resolved
wavelets. Their scales are a,...,b-1, their centres are p**a*r for
r=0,...,p**(c-a)-1, and each centre has p-1 frequencies. The labels are
``None`` for the constant mode and ``(c,r,j)`` for wavelets. The count is

.. math::

   1+\sum_{c=a}^{b-1}(p-1)p^{c-a}=p^{b-a}=N.

Disjoint supports give orthogonality. For nested supports at different
scales, the coarser wavelet is constant on the finer support, whose wavelet
has mean zero. At equal scale, the p-th-root character sum gives
orthogonality of distinct frequencies. Thus these N nonzero orthogonal
vectors span the N-dimensional declared function space. Decomposition
divides each inner product by its stored squared norm; synthesis and
weighted Parseval are exact.

.. doctest::

   >>> basis = KozyrevBasis(K, constancy_exponent=2)
   >>> len(basis.modes)
   9
   >>> f = PadicTestFunction(K, range(9), constancy_exponent=2)
   >>> coefficients = basis.decompose(f)
   >>> basis.reconstruct(coefficients) == f
   True
   >>> basis.energy(coefficients) == f.norm_squared()
   True

The whole-Qp Vladimirov derivative
-----------------------------------

For a positive integer order alpha the implementation uses

.. math::

   D^\alpha f=\mathcal F^{-1}(|\xi|_p^\alpha\widehat f),\qquad
   C_\alpha=\frac{p^\alpha-1}{1-p^{-1-\alpha}},\qquad
   D^\alpha f(x)=C_\alpha\int_{\mathbb Q_p}
      \frac{f(x)-f(y)}{|x-y|_p^{1+\alpha}}\,dy.

Kozyrev's paper uses the opposite Fourier sign; the radial multiplier is
even, so this operator agrees. His eigenvalue formula gives, in our scale,
``D^alpha(w)=p**(alpha*(c+1))*w``. See the primary source in :doc:`references`.
The native API currently requires positive integer alpha so all spectral
and kernel constants remain rational and exact.

For a general test function supported in B=p**a*Zp, a nonzero mass produces
an exterior tail. ``vladimirov`` returns a ``VladimirovImage`` containing
both the finite interior and this full analytic tail. For x outside B,
ultrametricity makes every distance to y in B equal to ``|x|``, hence

.. math::

   D^\alpha f(x)=-C_\alpha m\,|x|_p^{-1-\alpha},\qquad m=\int f(y)\,dy.

``image(point)`` evaluates either part. ``integral()`` and ``norm_squared()``
sum all exterior shells analytically using a convergent geometric series.
``as_test_function()`` succeeds exactly when the tail is zero. The full
derivative has integral zero, including when its interior alone does not.

.. doctest::

   >>> derivative = vladimirov(wave, order=2)
   >>> derivative.as_test_function() == wave*81
   True
   >>> indicator = PadicTestFunction.ball(K)
   >>> full = vladimirov(indicator)
   >>> full.interior(0).as_fraction(), full(Fraction(1, 3)).as_fraction()
   (Fraction(3, 4), Fraction(-1, 4))
   >>> full.integral().as_fraction()
   Fraction(0, 1)
   >>> full.as_test_function()
   Traceback (most recent call last):
       ...
   padic_tools.core.DomainError: whole-Qp Vladimirov image has a noncompact exterior tail

Conservative regional dynamics in a ball
-----------------------------------------

``BallVladimirov`` implements the regional integral over B. This choice
keeps jumps inside that ball and conserves its mass:

.. math::

   L_B f(x)=C_\alpha\int_B\frac{f(x)-f(y)}{|x-y|_p^{1+\alpha}}\,dy,
   \qquad x\in B.

For the zero extension of f, its whole-Qp derivative inside B is
``L_B(f)+kappa_a*f``. The escape rate follows by summing shells outside B:

.. math::

   \kappa_a=C_\alpha(1-p^{-1})\sum_{n=-\infty}^{a-1}p^{\alpha n}
      =\frac{(1-p^{-1})p^{\alpha a}}{1-p^{-1-\alpha}}.

This makes the boundary choice explicit. Killing jumps that leave B would
retain the escape term and lose interior mass; the regional operator
removes that term. It is the compact-ball integral discussed by Bikulov
and Zubarev, with our positive sign and factor C_alpha stated explicitly.

At frequency xi=p**(-b)*k, nonzero dual cosets have ``|xi|>=p**(a+1)``.
The whole-space Fourier multiplier is constant there. The regional
eigenvalues are zero for the constant mode and ``|xi|**alpha-kappa_a``
for all other modes. Resolved Kozyrev wavelets have eigenvalue
``p**(alpha*(c+1))-kappa_a``. Both forms are implemented and tested.

.. doctest::

   >>> operator = BallVladimirov(K, order=1, constancy_exponent=2)
   >>> operator.escape_rate
   Fraction(3, 4)
   >>> sorted(set(operator.eigenvalues))
   [Fraction(0, 1), Fraction(9, 4), Fraction(33, 4)]
   >>> operator.apply(f).integral().as_fraction()
   Fraction(0, 1)

The independent kernel matrix uses source-cell Haar weight p**(-b):

.. math::

   L_{jk}=-C_\alpha p^{-b}|p^a(j-k)|_p^{-1-\alpha}\quad(j\ne k),
   \qquad L_{jj}=-\sum_{k\ne j}L_{jk}.

It is symmetric, its off-diagonals are negative, and its rows sum to zero.
Therefore -L is a continuous-time Markov generator. Its quadratic form is
one half of the weighted sum of squared differences, so L is nonnegative.
These facts imply positivity, mass conservation and L2 contraction of
exp(-t*L). They are properties of this regional model with its declared
boundary, rather than assumptions about a sampled whole-Qp process.

Numerical diffusion and a reproducible example
----------------------------------------------

``operator.evolve(f, time, dps=...)`` solves ``du/dt=-L_B*u`` by multiplying
each finite Fourier mode by exp(-time*eigenvalue). Time must be exact
nonnegative rational data. The result is a ``DiffusionState`` whose values,
integral and squared norm are complex numerical approximations. Working
dps uses an isolated mpmath context and does not certify a complex error
ball. Exact operator data and initial mass remain available separately.

.. doctest::

   >>> initial = PadicTestFunction.ball(K, exponent=2)*9
   >>> initial = initial.on_lattice(support_exponent=0, constancy_exponent=2)
   >>> state = operator.evolve(initial, Fraction(1, 10), dps=60)
   >>> abs(state.integral()-1) < 10**-50
   True
   >>> min(value.real for value in state.values) > -10**-50
   True
   >>> state.norm_squared() < initial.norm_squared().complex(dps=60).real
   True

Run the complete comparison after installing the complex extra:

.. code-block:: console

   python examples/padic_diffusion.py

It uses p=3, order one and the same initial density ``9*1_(9*Z3)`` in three
configurations. The pair (a,b)=(0,2),(0,3) refines the grid inside the same
ball. The function already belongs to the coarser space, which is invariant
under L_B, so the two evolutions agree to numerical working accuracy.
The pair (-1,2) enlarges the ball and changes the regional boundary; it
conserves the same total mass but has a different stationary density.
The JSON output reports spectra, exact row-sum checks, density extrema,
mass errors, squared norms, and the refinement discrepancy.

.. _diffusion-reference-run:

Reference run
~~~~~~~~~~~~~~

The centre density in a reference run is shown below, rounded to six
decimals. Both cases start from exactly the same function. The unit-ball
b=2 and b=3 results agree to better than 10**(-50) in the checked run;
all reported mass errors are also below that threshold.

.. list-table:: Density at x=0
   :header-rows: 1

   * - Time
     - Ball Z3
     - Ball 3**(-1)*Z3
   * - 0
     - 9.000000
     - 9.000000
   * - 1/10
     - 5.226442
     - 4.972145
   * - 1
     - 1.212366
     - 0.777051
   * - 5
     - 1.000026
     - 0.349014

.. list-table:: Haar L2 norm squared
   :header-rows: 1

   * - Time
     - Ball Z3
     - Ball 3**(-1)*Z3
   * - 0
     - 9.000000
     - 9.000000
   * - 1/10
     - 3.427556
     - 3.103682
   * - 1
     - 1.022218
     - 0.490260
   * - 5
     - 1.000000
     - 0.333702

``on_lattice`` only enlarges the support container and/or refines constancy.
It represents exactly the same initial function. This refinement check is
not a universal convergence bound for arbitrary input or whole-space
diffusion. On a ball of volume p**(-a), the stationary density of unit mass
is p**a; the spectral gap determines its convergence rate.

Scope and validation
----------------------

Native wavelet bases and operators default to at most 256 cosets; explicit
``max_size`` can raise that bound, subject to the harmonic coefficient
backend's root-order limit. Bases occupy quadratic storage. Fourier,
matrix and synthesis operations are intended for small exact models.
General real-order operators, arbitrary initial-function approximation
bounds, multidimensional bases and whole-space numerical heat evolution
are subsequent extensions.

Validation includes exact Gram matrices, reconstruction, weighted Parseval,
kernel/Fourier agreement, wavelet eigenvalues and the full derivative's
exterior integral and norm. PARI/GP independently checks regional kernel
actions. Numerical heat solutions are checked against mpmath matrix
exponentials, positivity, mass, L2 contraction and refinement consistency.
