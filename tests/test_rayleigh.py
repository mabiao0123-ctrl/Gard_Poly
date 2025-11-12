"""Tests for Rayleigh polynomial functions."""

import numpy as np
import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from gard_poly.rayleigh import is_rayleigh, check_rayleigh


def test_is_rayleigh_basic():
    """Test basic Rayleigh polynomial detection."""
    # x^2 - 3x + 2 = (x-1)(x-2)
    # Roots: 1, 2
    # Derivative: 2x - 3, root at 1.5
    # Interlacing: 1 <= 1.5 <= 2 ✓
    assert is_rayleigh([1, -3, 2]), "x^2 - 3x + 2 should be Rayleigh"
    
    # x^2 - 1 = (x-1)(x+1)
    # Roots: -1, 1
    # Derivative: 2x, root at 0
    # Interlacing: -1 <= 0 <= 1 ✓
    assert is_rayleigh([1, 0, -1]), "x^2 - 1 should be Rayleigh"
    
    # Linear polynomial
    assert is_rayleigh([1, -1]), "Linear polynomials should be Rayleigh"
    
    # Constant polynomial
    assert is_rayleigh([5]), "Constant polynomial should be Rayleigh"
    
    print("✓ test_is_rayleigh_basic passed")


def test_is_rayleigh_non_real_roots():
    """Test that polynomials with non-real roots are not Rayleigh."""
    # x^2 + 1 = (x-i)(x+i)
    assert not is_rayleigh([1, 0, 1]), "x^2 + 1 should not be Rayleigh"
    
    print("✓ test_is_rayleigh_non_real_roots passed")


def test_is_rayleigh_higher_degree():
    """Test higher degree Rayleigh polynomials."""
    # (x-1)^2 = x^2 - 2x + 1
    # Roots: 1 (double root)
    # Derivative: 2x - 2, root at 1
    assert is_rayleigh([1, -2, 1]), "(x-1)^2 should be Rayleigh"
    
    # x^3 - 3x^2 + 2x = x(x-1)(x-2)
    # Roots: 0, 1, 2
    # Derivative: 3x^2 - 6x + 2
    # This has real roots that should interlace
    assert is_rayleigh([1, -3, 2, 0]), "x^3 - 3x^2 + 2x should be Rayleigh"
    
    print("✓ test_is_rayleigh_higher_degree passed")


def test_check_rayleigh():
    """Test detailed Rayleigh polynomial checking."""
    # Test Rayleigh polynomial
    result = check_rayleigh([1, -3, 2])
    assert result['is_rayleigh'], "x^2 - 3x + 2 should be Rayleigh"
    assert result['all_roots_real'], "All roots should be real"
    assert result['interlacing'], "Should have interlacing property"
    
    # Test non-Rayleigh polynomial
    result = check_rayleigh([1, 0, 1])
    assert not result['is_rayleigh'], "x^2 + 1 should not be Rayleigh"
    assert not result['all_roots_real'], "Not all roots are real"
    
    print("✓ test_check_rayleigh passed")


def test_check_rayleigh_verbose():
    """Test verbose output."""
    print("\nTesting verbose output for x^2 - 3x + 2:")
    result = check_rayleigh([1, -3, 2], verbose=True)
    assert result['is_rayleigh'], "Should be Rayleigh"
    
    print("\nTesting verbose output for x^2 + 1:")
    result = check_rayleigh([1, 0, 1], verbose=True)
    assert not result['is_rayleigh'], "Should not be Rayleigh"
    
    print("✓ test_check_rayleigh_verbose passed")


if __name__ == "__main__":
    test_is_rayleigh_basic()
    test_is_rayleigh_non_real_roots()
    test_is_rayleigh_higher_degree()
    test_check_rayleigh()
    test_check_rayleigh_verbose()
    print("\n✓ All Rayleigh tests passed!")
