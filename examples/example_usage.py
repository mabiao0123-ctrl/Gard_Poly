"""
Example usage of the gard_poly package.

This script demonstrates how to use the tools for testing
real stable, Rayleigh, and Lorentzian polynomials.
"""

import numpy as np
import sys
import os

# Add parent directory to path for running as script
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from gard_poly import (
    is_real_stable, check_real_stable,
    is_rayleigh, check_rayleigh,
    is_lorentzian, check_lorentzian
)


def main():
    """Run examples for all polynomial types."""
    
    print("=" * 70)
    print("GARD_POLY PACKAGE EXAMPLES")
    print("=" * 70)
    
    # Example 1: Real Stable Polynomials
    print("\n" + "=" * 70)
    print("Example 1: Real Stable Polynomials")
    print("=" * 70)
    
    print("\n1.1 Testing x^2 - 1 (roots at ±1, on real axis):")
    poly1 = [1, 0, -1]
    print(f"Polynomial coefficients: {poly1}")
    print(f"Is real stable: {is_real_stable(poly1)}")
    result = check_real_stable(poly1, verbose=True)
    
    print("\n1.2 Testing x^2 + 1 (roots at ±i, in complex plane):")
    poly2 = [1, 0, 1]
    print(f"Polynomial coefficients: {poly2}")
    print(f"Is real stable: {is_real_stable(poly2)}")
    result = check_real_stable(poly2, verbose=True)
    
    print("\n1.3 Testing (x-1)(x-2) = x^2 - 3x + 2:")
    poly3 = [1, -3, 2]
    print(f"Polynomial coefficients: {poly3}")
    print(f"Is real stable: {is_real_stable(poly3)}")
    
    # Example 2: Rayleigh Polynomials
    print("\n" + "=" * 70)
    print("Example 2: Rayleigh Polynomials")
    print("=" * 70)
    
    print("\n2.1 Testing x^2 - 3x + 2 (roots: 1, 2; derivative root: 1.5):")
    poly4 = [1, -3, 2]
    print(f"Polynomial coefficients: {poly4}")
    print(f"Is Rayleigh: {is_rayleigh(poly4)}")
    result = check_rayleigh(poly4, verbose=True)
    
    print("\n2.2 Testing x^2 + 1 (non-real roots):")
    poly5 = [1, 0, 1]
    print(f"Polynomial coefficients: {poly5}")
    print(f"Is Rayleigh: {is_rayleigh(poly5)}")
    result = check_rayleigh(poly5, verbose=True)
    
    print("\n2.3 Testing (x-1)^2 = x^2 - 2x + 1:")
    poly6 = [1, -2, 1]
    print(f"Polynomial coefficients: {poly6}")
    print(f"Is Rayleigh: {is_rayleigh(poly6)}")
    
    # Example 3: Lorentzian Polynomials
    print("\n" + "=" * 70)
    print("Example 3: Lorentzian Polynomials")
    print("=" * 70)
    
    print("\n3.1 Testing (x+1)^2 = x^2 + 2x + 1:")
    poly7 = [1, 2, 1]
    print(f"Polynomial coefficients: {poly7}")
    print(f"Is Lorentzian: {is_lorentzian(poly7)}")
    result = check_lorentzian(poly7, verbose=True)
    
    print("\n3.2 Testing x^2 + 1 (non-real roots):")
    poly8 = [1, 0, 1]
    print(f"Polynomial coefficients: {poly8}")
    print(f"Is Lorentzian: {is_lorentzian(poly8)}")
    result = check_lorentzian(poly8, verbose=True)
    
    print("\n3.3 Testing (x-1)^3 = x^3 - 3x^2 + 3x - 1:")
    poly9 = [1, -3, 3, -1]
    print(f"Polynomial coefficients: {poly9}")
    print(f"Is Lorentzian: {is_lorentzian(poly9)}")
    
    # Summary comparison
    print("\n" + "=" * 70)
    print("Summary Comparison")
    print("=" * 70)
    
    test_polynomials = [
        ([1, 0, -1], "x^2 - 1"),
        ([1, 0, 1], "x^2 + 1"),
        ([1, -3, 2], "x^2 - 3x + 2"),
        ([1, 2, 1], "x^2 + 2x + 1"),
        ([1, -2, 1], "x^2 - 2x + 1"),
    ]
    
    print("\n{:<20} {:<15} {:<15} {:<15}".format(
        "Polynomial", "Real Stable", "Rayleigh", "Lorentzian"
    ))
    print("-" * 70)
    
    for coeffs, name in test_polynomials:
        rs = "Yes" if is_real_stable(coeffs) else "No"
        ray = "Yes" if is_rayleigh(coeffs) else "No"
        lor = "Yes" if is_lorentzian(coeffs) else "No"
        print("{:<20} {:<15} {:<15} {:<15}".format(name, rs, ray, lor))
    
    print("\n" + "=" * 70)
    print("Examples completed!")
    print("=" * 70)


if __name__ == "__main__":
    main()
