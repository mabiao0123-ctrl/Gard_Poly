"""
Integration test demonstrating all package features.
"""

import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from gard_poly import (
    is_real_stable, check_real_stable,
    is_rayleigh, check_rayleigh,
    is_lorentzian, check_lorentzian,
    polynomial_roots, evaluate_polynomial
)


def test_integration():
    """Test all features together."""
    print("Running integration tests...")
    
    # Test 1: Polynomial that is real stable, Rayleigh, and Lorentzian
    print("\nTest 1: (x-1)(x-2) = x^2 - 3x + 2")
    poly1 = [1, -3, 2]
    
    assert is_real_stable(poly1), "Should be real stable"
    assert is_rayleigh(poly1), "Should be Rayleigh"
    assert is_lorentzian(poly1), "Should be Lorentzian"
    
    roots = polynomial_roots(poly1)
    assert len(roots) == 2, "Should have 2 roots"
    
    # Evaluate at x=1 (should be 0)
    val = evaluate_polynomial(poly1, 1)
    assert abs(val) < 1e-10, "Should be zero at x=1"
    
    print("  ✓ All properties verified")
    
    # Test 2: Polynomial that is none of the above
    print("\nTest 2: x^2 + 1 (complex roots)")
    poly2 = [1, 0, 1]
    
    assert not is_real_stable(poly2), "Should not be real stable"
    assert not is_rayleigh(poly2), "Should not be Rayleigh"
    assert not is_lorentzian(poly2), "Should not be Lorentzian"
    
    print("  ✓ All properties verified")
    
    # Test 3: Comprehensive check functions
    print("\nTest 3: Detailed check functions")
    poly3 = [1, 0, -1]
    
    rs_result = check_real_stable(poly3)
    assert rs_result['is_real_stable'], "Should be real stable"
    assert len(rs_result['upper_half_plane_roots']) == 0, "No upper half-plane roots"
    
    ray_result = check_rayleigh(poly3)
    assert ray_result['is_rayleigh'], "Should be Rayleigh"
    assert ray_result['interlacing'], "Should have interlacing property"
    
    lor_result = check_lorentzian(poly3)
    assert lor_result['is_lorentzian'], "Should be Lorentzian"
    assert lor_result['all_roots_real'], "All roots should be real"
    
    print("  ✓ All detailed checks passed")
    
    # Test 4: Edge cases
    print("\nTest 4: Edge cases")
    
    # Constant polynomial
    const_poly = [5]
    assert is_real_stable(const_poly), "Constant should be real stable"
    assert is_rayleigh(const_poly), "Constant should be Rayleigh"
    assert is_lorentzian(const_poly), "Constant should be Lorentzian"
    
    # Linear polynomial
    linear_poly = [1, -1]
    assert is_real_stable(linear_poly), "Linear should be real stable"
    assert is_rayleigh(linear_poly), "Linear should be Rayleigh"
    assert is_lorentzian(linear_poly), "Linear should be Lorentzian"
    
    print("  ✓ Edge cases handled correctly")
    
    print("\n" + "=" * 70)
    print("✓ ALL INTEGRATION TESTS PASSED!")
    print("=" * 70)


if __name__ == "__main__":
    test_integration()
