#!/usr/bin/env python3
"""
Run all tests for the gard_poly package.
"""

import sys
import subprocess


def run_test_file(test_file):
    """Run a single test file and return True if successful."""
    print(f"\n{'=' * 70}")
    print(f"Running {test_file}")
    print('=' * 70)
    
    result = subprocess.run(
        [sys.executable, test_file],
        capture_output=False
    )
    
    return result.returncode == 0


def main():
    """Run all test files."""
    test_files = [
        'tests/test_utils.py',
        'tests/test_real_stable.py',
        'tests/test_rayleigh.py',
        'tests/test_lorentzian.py',
    ]
    
    print("=" * 70)
    print("GARD_POLY TEST SUITE")
    print("=" * 70)
    
    all_passed = True
    for test_file in test_files:
        if not run_test_file(test_file):
            all_passed = False
    
    print("\n" + "=" * 70)
    if all_passed:
        print("✓ ALL TESTS PASSED!")
    else:
        print("✗ SOME TESTS FAILED")
        sys.exit(1)
    print("=" * 70)


if __name__ == "__main__":
    main()
