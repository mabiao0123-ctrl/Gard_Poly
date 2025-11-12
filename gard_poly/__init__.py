"""
Gard_Poly: Tools for testing real stable, Rayleigh, and Lorentzian polynomials.
"""

from .real_stable import is_real_stable, check_real_stable
from .rayleigh import is_rayleigh, check_rayleigh
from .lorentzian import is_lorentzian, check_lorentzian
from .utils import polynomial_roots, evaluate_polynomial

__version__ = "0.1.0"

__all__ = [
    "is_real_stable",
    "check_real_stable",
    "is_rayleigh",
    "check_rayleigh",
    "is_lorentzian",
    "check_lorentzian",
    "polynomial_roots",
    "evaluate_polynomial",
]
