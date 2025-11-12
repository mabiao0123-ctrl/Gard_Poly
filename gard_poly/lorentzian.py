"""
Tools for testing Lorentzian polynomials.

A homogeneous polynomial is Lorentzian if it satisfies certain conditions
related to the Hessian matrix and hyperbolicity. For simplicity, we implement
basic checks for univariate cases.
"""

import numpy as np
from .utils import polynomial_roots, all_real


def is_lorentzian(coefficients, tolerance=1e-5):
    """
    Check if a univariate polynomial is Lorentzian.
    
    For univariate polynomials, a polynomial is Lorentzian if:
    1. All roots are real
    2. The polynomial satisfies log-concavity conditions on its coefficients
    
    For a polynomial with non-negative coefficients a_0, a_1, ..., a_n,
    log-concavity means: a_i^2 >= a_{i-1} * a_{i+1} for all i.
    
    Parameters:
    -----------
    coefficients : array-like
        Polynomial coefficients in descending order of powers.
    tolerance : float, optional
        Numerical tolerance for comparisons (default: 1e-10).
    
    Returns:
    --------
    bool
        True if the polynomial is Lorentzian (in the univariate sense), False otherwise.
    
    Examples:
    ---------
    >>> is_lorentzian([1, 2, 1])  # (x+1)^2 = x^2 + 2x + 1
    True
    >>> is_lorentzian([1, -2, 1])  # (x-1)^2 = x^2 - 2x + 1
    True
    """
    coefficients = np.array(coefficients, dtype=np.float64)
    
    # Remove leading zeros
    while len(coefficients) > 1 and np.abs(coefficients[0]) < tolerance:
        coefficients = coefficients[1:]
    
    if len(coefficients) <= 2:
        # Linear and constant polynomials are trivially Lorentzian
        return True
    
    # Get roots
    roots = polynomial_roots(coefficients, tolerance)
    
    # Check if all roots are real
    if len(roots) > 0 and not all_real(roots, tolerance):
        return False
    
    # For polynomials with all real roots, we also check a form of log-concavity
    # This is a simplified check - true Lorentzian polynomials in multiple
    # variables have more complex requirements
    
    # Check log-concavity of absolute values of coefficients
    abs_coeffs = np.abs(coefficients)
    
    # Skip if any coefficient is too close to zero (to avoid division issues)
    if np.any(abs_coeffs < tolerance):
        # If we have zeros in coefficients, just rely on real roots check
        return True
    
    # Check a_i^2 >= a_{i-1} * a_{i+1} for middle coefficients
    for i in range(1, len(abs_coeffs) - 1):
        if abs_coeffs[i] ** 2 < abs_coeffs[i-1] * abs_coeffs[i+1] - tolerance:
            return False
    
    return True


def check_lorentzian(coefficients, tolerance=1e-5, verbose=False):
    """
    Check if a polynomial is Lorentzian and provide detailed information.
    
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
        - 'is_lorentzian': bool, whether the polynomial is Lorentzian
        - 'roots': numpy.ndarray, the roots of the polynomial
        - 'all_roots_real': bool, whether all roots are real
        - 'log_concave': bool, whether coefficients satisfy log-concavity
    
    Examples:
    ---------
    >>> result = check_lorentzian([1, 2, 1])
    >>> result['is_lorentzian']
    True
    """
    coefficients = np.array(coefficients, dtype=np.float64)
    
    # Remove leading zeros
    while len(coefficients) > 1 and np.abs(coefficients[0]) < tolerance:
        coefficients = coefficients[1:]
    
    # Get roots
    roots = polynomial_roots(coefficients, tolerance)
    all_roots_real = all_real(roots, tolerance) if len(roots) > 0 else True
    
    # Check log-concavity
    abs_coeffs = np.abs(coefficients)
    log_concave = True
    
    if len(coefficients) > 2:
        for i in range(1, len(abs_coeffs) - 1):
            if abs_coeffs[i] < tolerance:
                continue  # Skip zero coefficients
            if abs_coeffs[i] ** 2 < abs_coeffs[i-1] * abs_coeffs[i+1] - tolerance:
                log_concave = False
                break
    
    is_lorentzian_poly = all_roots_real and log_concave
    
    result = {
        'is_lorentzian': is_lorentzian_poly,
        'roots': roots,
        'all_roots_real': all_roots_real,
        'log_concave': log_concave,
    }
    
    if verbose:
        print(f"Lorentzian Polynomial Check:")
        print(f"  Coefficients: {coefficients}")
        print(f"  Roots: {roots}")
        print(f"  All roots real: {all_roots_real}")
        print(f"  Log-concave coefficients: {log_concave}")
        print(f"  Is Lorentzian: {is_lorentzian_poly}")
    
    return result
