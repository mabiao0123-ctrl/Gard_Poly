"""
Utility functions for polynomial operations.
"""

import numpy as np


def polynomial_roots(coefficients, tolerance=1e-5):
    """
    Compute roots of a polynomial given its coefficients.
    
    Parameters:
    -----------
    coefficients : array-like
        Polynomial coefficients in descending order of powers.
        For example, [1, -3, 2] represents x^2 - 3x + 2.
    tolerance : float, optional
        Numerical tolerance for comparisons (default: 1e-10).
    
    Returns:
    --------
    numpy.ndarray
        Array of polynomial roots (may be complex).
    
    Examples:
    ---------
    >>> roots = polynomial_roots([1, -3, 2])  # x^2 - 3x + 2 = (x-1)(x-2)
    >>> np.sort(roots)
    array([1., 2.])
    """
    coefficients = np.array(coefficients, dtype=np.complex128)
    
    # Remove leading zeros
    while len(coefficients) > 1 and np.abs(coefficients[0]) < tolerance:
        coefficients = coefficients[1:]
    
    if len(coefficients) == 0:
        return np.array([])
    
    if len(coefficients) == 1:
        if np.abs(coefficients[0]) < tolerance:
            return np.array([])  # Zero polynomial
        else:
            return np.array([])  # Constant polynomial has no roots
    
    return np.roots(coefficients)


def evaluate_polynomial(coefficients, x):
    """
    Evaluate polynomial at given value(s).
    
    Parameters:
    -----------
    coefficients : array-like
        Polynomial coefficients in descending order of powers.
    x : complex or array-like
        Value(s) at which to evaluate the polynomial.
    
    Returns:
    --------
    complex or numpy.ndarray
        Polynomial value(s) at x.
    
    Examples:
    ---------
    >>> evaluate_polynomial([1, -3, 2], 0)  # x^2 - 3x + 2 at x=0
    (2+0j)
    >>> evaluate_polynomial([1, -3, 2], 1)  # x^2 - 3x + 2 at x=1
    0j
    """
    coefficients = np.array(coefficients, dtype=np.complex128)
    x = np.asarray(x, dtype=np.complex128)
    return np.polyval(coefficients, x)


def is_real(z, tolerance=1e-5):
    """
    Check if a complex number is effectively real within tolerance.
    
    Parameters:
    -----------
    z : complex or array-like
        Number(s) to check.
    tolerance : float, optional
        Numerical tolerance (default: 1e-10).
    
    Returns:
    --------
    bool or numpy.ndarray
        True if imaginary part is within tolerance of zero.
    """
    z = np.asarray(z, dtype=np.complex128)
    return np.abs(z.imag) < tolerance


def all_real(roots, tolerance=1e-5):
    """
    Check if all roots are real within tolerance.
    
    Parameters:
    -----------
    roots : array-like
        Array of complex numbers.
    tolerance : float, optional
        Numerical tolerance (default: 1e-10).
    
    Returns:
    --------
    bool
        True if all roots are real within tolerance.
    """
    if len(roots) == 0:
        return True
    roots = np.asarray(roots, dtype=np.complex128)
    return np.all(is_real(roots, tolerance))
