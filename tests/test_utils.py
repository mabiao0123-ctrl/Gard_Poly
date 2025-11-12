"""Tests for utility functions."""

import numpy as np
import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from gard_poly.utils import (
    polynomial_roots, 
    evaluate_polynomial, 
    is_real, 
    all_real
)


def test_polynomial_roots():
    """Test polynomial root finding."""
    # Test x^2 - 3x + 2 = (x-1)(x-2)
    roots = polynomial_roots([1, -3, 2])
    roots_sorted = np.sort(roots.real)
    assert np.allclose(roots_sorted, [1, 2]), f"Expected [1, 2], got {roots_sorted}"
    
    # Test x^2 + 1 = (x-i)(x+i)
    roots = polynomial_roots([1, 0, 1])
    assert len(roots) == 2
    assert np.allclose(np.sort(roots.imag), [-1, 1]), "Expected roots at ±i"
    
    # Test constant polynomial
    roots = polynomial_roots([5])
    assert len(roots) == 0, "Constant polynomial should have no roots"
    
    print("✓ test_polynomial_roots passed")


def test_evaluate_polynomial():
    """Test polynomial evaluation."""
    # Test x^2 - 3x + 2 at x=0
    val = evaluate_polynomial([1, -3, 2], 0)
    assert np.allclose(val, 2), f"Expected 2, got {val}"
    
    # Test x^2 - 3x + 2 at x=1
    val = evaluate_polynomial([1, -3, 2], 1)
    assert np.allclose(val, 0), f"Expected 0, got {val}"
    
    # Test x^2 - 3x + 2 at x=2
    val = evaluate_polynomial([1, -3, 2], 2)
    assert np.allclose(val, 0), f"Expected 0, got {val}"
    
    print("✓ test_evaluate_polynomial passed")


def test_is_real():
    """Test real number checking."""
    assert is_real(1.0 + 0j), "1.0 should be real"
    assert is_real(5.0), "5.0 should be real"
    assert not is_real(1.0 + 1j), "1.0 + 1j should not be real"
    
    # Test with tolerance
    assert is_real(1.0 + 1e-12j, tolerance=1e-10), "Should be real within tolerance"
    
    print("✓ test_is_real passed")


def test_all_real():
    """Test checking if all roots are real."""
    assert all_real([1.0, 2.0, 3.0]), "All should be real"
    assert all_real([1.0 + 1e-12j, 2.0]), "Should be real within tolerance"
    assert not all_real([1.0, 2.0j]), "2j is not real"
    assert all_real([]), "Empty list should be considered all real"
    
    print("✓ test_all_real passed")


if __name__ == "__main__":
    test_polynomial_roots()
    test_evaluate_polynomial()
    test_is_real()
    test_all_real()
    print("\n✓ All utility tests passed!")
