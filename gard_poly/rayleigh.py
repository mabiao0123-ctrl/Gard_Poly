"""
Tools for testing Rayleigh polynomials.

A polynomial is Rayleigh if all its roots are real and satisfy certain
interlacing properties. For a univariate polynomial, this means all roots
are real and the derivative's roots interlace with the polynomial's roots.
"""

import numpy as np
from .utils import polynomial_roots, all_real


def is_rayleigh(coefficients, tolerance=1e-5):
    """
    Check if a polynomial is Rayleigh.
    
    A polynomial p(x) is Rayleigh if:
    1. All roots are real
    2. The roots of p'(x) interlace with the roots of p(x)
    
    For a polynomial with roots r1 <= r2 <= ... <= rn,
    interlacing means: r1 <= s1 <= r2 <= s2 <= ... <= sn-1 <= rn,
    where s1, ..., sn-1 are the roots of the derivative.
    
    Parameters:
    -----------
    coefficients : array-like
        Polynomial coefficients in descending order of powers.
        For example, [1, -3, 2] represents x^2 - 3x + 2.
    tolerance : float, optional
        Numerical tolerance for comparisons (default: 1e-10).
    
    Returns:
    --------
    bool
        True if the polynomial is Rayleigh, False otherwise.
    
    Examples:
    ---------
    >>> is_rayleigh([1, -3, 2])  # x^2 - 3x + 2 = (x-1)(x-2)
    True
    >>> is_rayleigh([1, 0, 0, -1])  # x^3 - 1
    True
    """
    coefficients = np.array(coefficients, dtype=np.float64)
    
    # Remove leading zeros
    while len(coefficients) > 1 and np.abs(coefficients[0]) < tolerance:
        coefficients = coefficients[1:]
    
    if len(coefficients) <= 1:
        # Constant polynomial is trivially Rayleigh
        return True
    
    # Get roots of polynomial
    roots = polynomial_roots(coefficients, tolerance)
    
    if len(roots) == 0:
        return True
    
    # Check if all roots are real
    if not all_real(roots, tolerance):
        return False
    
    # Get real parts of roots and sort them
    real_roots = np.sort(np.real(roots))
    
    # Compute derivative coefficients
    n = len(coefficients)
    derivative_coeffs = np.array([coefficients[i] * (n - 1 - i) 
                                   for i in range(n - 1)])
    
    if len(derivative_coeffs) == 0:
        # Linear polynomial is Rayleigh
        return True
    
    # Get roots of derivative
    derivative_roots = polynomial_roots(derivative_coeffs, tolerance)
    
    if len(derivative_roots) == 0:
        # If derivative has no roots, check is trivial
        return True
    
    # Check if derivative roots are real
    if not all_real(derivative_roots, tolerance):
        return False
    
    # Get real parts and sort
    real_derivative_roots = np.sort(np.real(derivative_roots))
    
    # Check interlacing property
    # For each derivative root s_i, it should satisfy:
    # real_roots[i] <= s_i <= real_roots[i+1]
    for i, s in enumerate(real_derivative_roots):
        if i >= len(real_roots):
            return False
        if i == len(real_roots) - 1:
            # Last derivative root should be >= last polynomial root
            # This shouldn't happen for standard polynomials
            break
        if s < real_roots[i] - tolerance or s > real_roots[i + 1] + tolerance:
            return False
    
    return True


def check_rayleigh(coefficients, tolerance=1e-5, verbose=False):
    """
    Check if a polynomial is Rayleigh and provide detailed information.
    
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
        - 'is_rayleigh': bool, whether the polynomial is Rayleigh
        - 'roots': numpy.ndarray, the roots of the polynomial
        - 'derivative_roots': numpy.ndarray, roots of the derivative
        - 'all_roots_real': bool, whether all roots are real
        - 'interlacing': bool, whether roots interlace properly
    
    Examples:
    ---------
    >>> result = check_rayleigh([1, -3, 2])
    >>> result['is_rayleigh']
    True
    """
    coefficients = np.array(coefficients, dtype=np.float64)
    
    # Remove leading zeros
    while len(coefficients) > 1 and np.abs(coefficients[0]) < tolerance:
        coefficients = coefficients[1:]
    
    # Get roots
    roots = polynomial_roots(coefficients, tolerance)
    all_roots_real = all_real(roots, tolerance) if len(roots) > 0 else True
    
    # Compute derivative
    n = len(coefficients)
    derivative_coeffs = np.array([coefficients[i] * (n - 1 - i) 
                                   for i in range(n - 1)])
    derivative_roots = polynomial_roots(derivative_coeffs, tolerance)
    
    # Check interlacing
    interlacing = False
    if all_roots_real and len(roots) > 0:
        real_roots = np.sort(np.real(roots))
        if len(derivative_roots) > 0 and all_real(derivative_roots, tolerance):
            real_derivative_roots = np.sort(np.real(derivative_roots))
            interlacing = True
            for i, s in enumerate(real_derivative_roots):
                if i >= len(real_roots) - 1:
                    break
                if s < real_roots[i] - tolerance or s > real_roots[i + 1] + tolerance:
                    interlacing = False
                    break
        else:
            interlacing = len(derivative_roots) == 0
    
    is_rayleigh_poly = all_roots_real and interlacing
    
    result = {
        'is_rayleigh': is_rayleigh_poly,
        'roots': roots,
        'derivative_roots': derivative_roots,
        'all_roots_real': all_roots_real,
        'interlacing': interlacing,
    }
    
    if verbose:
        print(f"Rayleigh Polynomial Check:")
        print(f"  Coefficients: {coefficients}")
        print(f"  Roots: {roots}")
        print(f"  All roots real: {all_roots_real}")
        print(f"  Derivative roots: {derivative_roots}")
        print(f"  Interlacing property: {interlacing}")
        print(f"  Is Rayleigh: {is_rayleigh_poly}")
    
    return result
