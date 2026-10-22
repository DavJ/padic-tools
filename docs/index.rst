padic-tools
===========

Numerical p-adic analysis, exact symbolic algebra, and explicit algebraic
bridges to complex numbers, quaternions, biquaternions and geometric algebra.
Version 0.1 keeps the coefficient data, analytic domain and approximation
precision visible in the API.

Install
-------

.. code-block:: console

   python -m pip install .
   python -m pip install '.[symbolic,complex,docs]'

The numerical core has no external dependency. SymPy, mpmath and Sphinx are
optional. HTML documentation is built with the ``docs`` extra; the complete
API reference includes the optional SymPy adapter and therefore also needs
the ``symbolic`` extra, or the combined ``dev`` extra.

First calculation
-----------------

.. doctest::

   >>> K = Qp(5, precision=20)
   >>> x = K(Fraction(1, 3))
   >>> (3*x).is_congruent(1)
   True
   >>> log(exp(K(5))).is_congruent(5)
   True
   >>> zeta_p(-3, K, branch=0).is_congruent(Fraction(-31, 30))
   True

Read :doc:`bridges` for the included complex/p-adic construction and its
Clifford interpretation. Read :doc:`functions` before comparing zeta values
across packages: the prime, branch and Euler factor must match.

.. toctree::
   :maxdepth: 1
   :caption: Guide

   Quick start <tutorial>
   Precision <precision>
   Special functions <functions>
   Integral transforms <transforms>
   Bounded measures <measures>
   Mellin of rational measures <mellin>
   Taylor expansions <taylor>
   Symbolic algebra <symbolic>
   Algebraic bridges <bridges>
   Architecture <architecture>
   Roadmap <roadmap>
   References <references>
   API reference <api>

.. toctree::
   :hidden:

   cesky
