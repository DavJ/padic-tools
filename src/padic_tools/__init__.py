"""p-adic numerical analysis and explicit exact algebraic bridges.

The native core has no third-party dependency. Optional symbolic, complex and
documentation extras are imported only when their corresponding APIs are used.
"""
from .core import DomainError, Padic, PadicField, PrecisionError, Qp, exact, valuation
from .analytic import cos, exp, log, log1p, power, sin, sqrt, teichmuller
from .special import bernoulli, theta3, theta4, zeta_euler_value, zeta_p, zeta_special
from .symbolic import Polynomial, PowerSeries, exp_series, hensel_lift, log1p_series, theta3_series
from .algebra import (
    Biquaternions, CliffordAlgebra, Gaussian, Hamilton, Multivector,
    QuadraticAlgebra, QuadraticElement, Quaternion, QuaternionAlgebra,
)
from .bridges import GaussianBridge, QuaternionBridge, RationalBridge, riemann_zeta

__version__ = "0.1.0"
__all__ = [
    "Qp", "PadicField", "Padic", "DomainError", "PrecisionError", "exact", "valuation",
    "exp", "log", "log1p", "sin", "cos", "sqrt", "power", "teichmuller",
    "bernoulli", "theta3", "theta4", "zeta_special", "zeta_euler_value", "zeta_p",
    "Polynomial", "PowerSeries", "hensel_lift", "exp_series", "log1p_series", "theta3_series",
    "QuadraticAlgebra", "QuadraticElement", "Gaussian", "QuaternionAlgebra", "Quaternion",
    "Hamilton", "Biquaternions", "CliffordAlgebra", "Multivector",
    "RationalBridge", "GaussianBridge", "QuaternionBridge", "riemann_zeta",
]
