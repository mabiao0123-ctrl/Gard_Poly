"""Tests for real stable polynomial functions."""

import numpy as np
import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from gard_poly.real_stable import is_real_stable, check_real_stable


def test_is_real_stable_basic():
    """Test basic real stable polynomial detection."""
    # x^2 - 1 = (x-1)(x+1) - roots at ±1 (on real axis)
    assert is_real_stable([1, 0, -1]), "x^2 - 1 should be real stable"
    
    # x^2 + 1 = (x-i)(x+i) - roots at ±i (in upper/lower half-plane)
    assert not is_real_stable([1, 0, 1]), "x^2 + 1 should not be real stable"
    
    # x - 1 - root at 1 (on real axis)
    assert is_real_stable([1, -1]), "x - 1 should be real stable"
    
    # Constant polynomial
    assert is_real_stable([5]), "Constant polynomial should be real stable"
    
    print("✓ test_is_real_stable_basic passed")


def test_is_real_stable_complex_coefficients():
    """Test that polynomials with complex coefficients are not real stable."""
    # Polynomial with complex coefficient
    assert not is_real_stable([1 + 1j, 0, -1]), "Complex coefficients should not be real stable"
    
    print("✓ test_is_real_stable_complex_coefficients passed")


def test_is_real_stable_multiple_roots():
    """Test polynomials with multiple real roots."""
    # (x-1)^2 = x^2 - 2x + 1
    assert is_real_stable([1, -2, 1]), "(x-1)^2 should be real stable"
    
    # (x-1)(x-2)(x-3) = x^3 - 6x^2 + 11x - 6
    assert is_real_stable([1, -6, 11, -6]), "(x-1)(x-2)(x-3) should be real stable"
    
    print("✓ test_is_real_stable_multiple_roots passed")


def test_check_real_stable():
    """Test detailed real stable polynomial checking."""
    # Test real stable polynomial
    result = check_real_stable([1, 0, -1])
    assert result['is_real_stable'], "x^2 - 1 should be real stable"
    assert result['real_coefficients'], "Coefficients should be real"
    assert len(result['upper_half_plane_roots']) == 0, "Should have no upper half-plane roots"
    
    # Test non-real stable polynomial
    result = check_real_stable([1, 0, 1])
    assert not result['is_real_stable'], "x^2 + 1 should not be real stable"
    assert len(result['upper_half_plane_roots']) > 0, "Should have upper half-plane roots"
    
    print("✓ test_check_real_stable passed")


def test_check_real_stable_verbose():
    """Test verbose output."""
    print("\nTesting verbose output for x^2 - 1:")
    result = check_real_stable([1, 0, -1], verbose=True)
    assert result['is_real_stable'], "Should be real stable"
    
    print("\nTesting verbose output for x^2 + 1:")
    result = check_real_stable([1, 0, 1], verbose=True)
    assert not result['is_real_stable'], "Should not be real stable"
    
    print("✓ test_check_real_stable_verbose passed")


if __name__ == "__main__":
    test_is_real_stable_basic()
    test_is_real_stable_complex_coefficients()
    test_is_real_stable_multiple_roots()
    test_check_real_stable()
    test_check_real_stable_verbose()
    print("\n✓ All real stable tests passed!")
