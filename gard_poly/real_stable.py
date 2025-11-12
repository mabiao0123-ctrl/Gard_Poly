"""
Tools for testing real stable polynomials.

A polynomial is real stable if it has no zeros in the upper half-plane
(i.e., all zeros have non-positive imaginary parts).
"""

import numpy as np
from .utils import polynomial_roots, evaluate_polynomial


def is_real_stable(coefficients, tolerance=1e-5):
    """
    Check if a polynomial is real stable.
    
    A polynomial p(z) is real stable if:
    1. All coefficients are real
    2. p(z) ≠ 0 for all z with Im(z) > 0
    
    This is equivalent to checking that all roots have non-positive imaginary parts.
    
    Parameters:
    -----------
    coefficients : array-like
        Polynomial coefficients in descending order of powers.
        For example, [1, 0, -1] represents z^2 - 1.
    tolerance : float, optional
        Numerical tolerance for comparisons (default: 1e-10).
    
    Returns:
    --------
    bool
        True if the polynomial is real stable, False otherwise.
    
    Examples:
    ---------
    >>> is_real_stable([1, 0, -1])  # z^2 - 1 = (z-1)(z+1)
    True
    >>> is_real_stable([1, 0, 1])   # z^2 + 1 = (z-i)(z+i)
    False
    """
    coefficients = np.array(coefficients, dtype=np.complex128)
    
    # Check if all coefficients are real
    if not np.all(np.abs(coefficients.imag) < tolerance):
        return False
    
    # Get roots
    roots = polynomial_roots(coefficients, tolerance)
    
    if len(roots) == 0:
        # Constant polynomial is real stable if it's real
        return True
    
    # Check if all roots have non-positive imaginary parts
    # A root is in upper half-plane if Im(z) > tolerance
    # We accept roots on real axis (|Im(z)| <= tolerance) or lower half-plane (Im(z) < -tolerance)
    for root in roots:
        if root.imag > tolerance:
            return False
    return True


def check_real_stable(coefficients, tolerance=1e-5, verbose=False):
    """
    Check if a polynomial is real stable and provide detailed information.
    
    Parameters:
    -----------
    coefficients : array-like
        Polynomial coefficients in descending order of powers.
    tolerance : float, optional
        Numerical tolerance for comparisons (default: 1e-10).
    verbose : bool, optional
        If True, print detailed information (default: False).
    
    Returns:
    --------
    dict
        Dictionary containing:
        - 'is_real_stable': bool, whether the polynomial is real stable
        - 'roots': numpy.ndarray, the roots of the polynomial
        - 'real_coefficients': bool, whether all coefficients are real
        - 'upper_half_plane_roots': list, roots in upper half-plane (if any)
    
    Examples:
    ---------
    >>> result = check_real_stable([1, 0, -1])
    >>> result['is_real_stable']
    True
    """
    coefficients = np.array(coefficients, dtype=np.complex128)
    
    # Check if coefficients are real
    real_coefficients = np.all(np.abs(coefficients.imag) < tolerance)
    
    # Get roots
    roots = polynomial_roots(coefficients, tolerance)
    
    # Find roots in upper half-plane
    upper_half_plane_roots = []
    if len(roots) > 0:
        upper_half_plane_roots = [r for r in roots if r.imag > tolerance]
    
    is_stable = real_coefficients and len(upper_half_plane_roots) == 0
    
    result = {
        'is_real_stable': is_stable,
        'roots': roots,
        'real_coefficients': real_coefficients,
        'upper_half_plane_roots': upper_half_plane_roots,
    }
    
    if verbose:
        print(f"Real Stable Polynomial Check:")
        print(f"  Coefficients: {coefficients}")
        print(f"  Real coefficients: {real_coefficients}")
        print(f"  Roots: {roots}")
        print(f"  Roots in upper half-plane: {upper_half_plane_roots}")
        print(f"  Is real stable: {is_stable}")
    
    return result
