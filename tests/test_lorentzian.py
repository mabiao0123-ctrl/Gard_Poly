"""Tests for Lorentzian polynomial functions."""

import numpy as np
import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from gard_poly.lorentzian import is_lorentzian, check_lorentzian


def test_is_lorentzian_basic():
    """Test basic Lorentzian polynomial detection."""
    # (x+1)^2 = x^2 + 2x + 1
    # All roots real, coefficients: [1, 2, 1]
    # Log-concavity: 2^2 = 4 >= 1*1 = 1 ✓
    assert is_lorentzian([1, 2, 1]), "(x+1)^2 should be Lorentzian"
    
    # (x-1)^2 = x^2 - 2x + 1
    # All roots real, absolute coefficients: [1, 2, 1]
    # Log-concavity: 2^2 = 4 >= 1*1 = 1 ✓
    assert is_lorentzian([1, -2, 1]), "(x-1)^2 should be Lorentzian"
    
    # Linear polynomial
    assert is_lorentzian([1, -1]), "Linear polynomials should be Lorentzian"
    
    # Constant polynomial
    assert is_lorentzian([5]), "Constant polynomial should be Lorentzian"
    
    print("✓ test_is_lorentzian_basic passed")


def test_is_lorentzian_non_real_roots():
    """Test that polynomials with non-real roots are not Lorentzian."""
    # x^2 + 1 = (x-i)(x+i)
    assert not is_lorentzian([1, 0, 1]), "x^2 + 1 should not be Lorentzian"
    
    print("✓ test_is_lorentzian_non_real_roots passed")


def test_is_lorentzian_real_roots():
    """Test polynomials with real roots."""
    # x^2 - 1 = (x-1)(x+1)
    # Roots: -1, 1 (real)
    # Coefficients: [1, 0, -1], absolute values: [1, 0, 1]
    # With zero coefficient, should still pass
    assert is_lorentzian([1, 0, -1]), "x^2 - 1 should be Lorentzian"
    
    # x^2 - 3x + 2 = (x-1)(x-2)
    # Roots: 1, 2 (real)
    # Coefficients: [1, -3, 2], absolute values: [1, 3, 2]
    # Log-concavity: 3^2 = 9 >= 1*2 = 2 ✓
    assert is_lorentzian([1, -3, 2]), "x^2 - 3x + 2 should be Lorentzian"
    
    print("✓ test_is_lorentzian_real_roots passed")


def test_check_lorentzian():
    """Test detailed Lorentzian polynomial checking."""
    # Test Lorentzian polynomial
    result = check_lorentzian([1, 2, 1])
    assert result['is_lorentzian'], "(x+1)^2 should be Lorentzian"
    assert result['all_roots_real'], "All roots should be real"
    assert result['log_concave'], "Should be log-concave"
    
    # Test non-Lorentzian polynomial
    result = check_lorentzian([1, 0, 1])
    assert not result['is_lorentzian'], "x^2 + 1 should not be Lorentzian"
    assert not result['all_roots_real'], "Not all roots are real"
    
    print("✓ test_check_lorentzian passed")


def test_check_lorentzian_verbose():
    """Test verbose output."""
    print("\nTesting verbose output for (x+1)^2:")
    result = check_lorentzian([1, 2, 1], verbose=True)
    assert result['is_lorentzian'], "Should be Lorentzian"
    
    print("\nTesting verbose output for x^2 + 1:")
    result = check_lorentzian([1, 0, 1], verbose=True)
    assert not result['is_lorentzian'], "Should not be Lorentzian"
    
    print("✓ test_check_lorentzian_verbose passed")


def test_is_lorentzian_higher_degree():
    """Test higher degree polynomials."""
    # x^3 - 6x^2 + 11x - 6 = (x-1)(x-2)(x-3)
    # All simple real roots, easier for numerical computation
    assert is_lorentzian([1, -6, 11, -6]), "(x-1)(x-2)(x-3) should be Lorentzian"
    
    print("✓ test_is_lorentzian_higher_degree passed")


if __name__ == "__main__":
    test_is_lorentzian_basic()
    test_is_lorentzian_non_real_roots()
    test_is_lorentzian_real_roots()
    test_check_lorentzian()
    test_check_lorentzian_verbose()
    test_is_lorentzian_higher_degree()
    print("\n✓ All Lorentzian tests passed!")
