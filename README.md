# Gard_Poly

Package for testing and managing polynomials such as real stable, Lorentzian, Rayleigh, Garding etc.

## Features

This package provides tools to test various types of special polynomials:

- **Real Stable Polynomials**: Polynomials with no zeros in the upper half-plane
- **Rayleigh Polynomials**: Polynomials with all real roots where the derivative's roots interlace with the polynomial's roots
- **Lorentzian Polynomials**: Polynomials satisfying log-concavity conditions (univariate case)

## Installation

```bash
pip install numpy scipy
```

Then install the package:
```bash
pip install -e .
```

Or install from source:
```bash
python setup.py install
```

## Usage

### Real Stable Polynomials

```python
from gard_poly import is_real_stable, check_real_stable

# Test if x^2 - 1 is real stable (roots at ±1)
coeffs = [1, 0, -1]
print(is_real_stable(coeffs))  # True

# Get detailed information
result = check_real_stable(coeffs, verbose=True)
print(result)
```

### Rayleigh Polynomials

```python
from gard_poly import is_rayleigh, check_rayleigh

# Test if x^2 - 3x + 2 is Rayleigh (roots: 1, 2; derivative root: 1.5)
coeffs = [1, -3, 2]
print(is_rayleigh(coeffs))  # True

# Get detailed information
result = check_rayleigh(coeffs, verbose=True)
print(result)
```

### Lorentzian Polynomials

```python
from gard_poly import is_lorentzian, check_lorentzian

# Test if (x+1)^2 is Lorentzian
coeffs = [1, 2, 1]
print(is_lorentzian(coeffs))  # True

# Get detailed information
result = check_lorentzian(coeffs, verbose=True)
print(result)
```

## Examples

Run the example script to see demonstrations:

```bash
python examples/example_usage.py
```

## Testing

Run the test suite:

```bash
python tests/test_utils.py
python tests/test_real_stable.py
python tests/test_rayleigh.py
python tests/test_lorentzian.py
```

## API Reference

### Real Stable Functions

- `is_real_stable(coefficients, tolerance=1e-5)`: Check if a polynomial is real stable
- `check_real_stable(coefficients, tolerance=1e-5, verbose=False)`: Get detailed real stable check

### Rayleigh Functions

- `is_rayleigh(coefficients, tolerance=1e-5)`: Check if a polynomial is Rayleigh
- `check_rayleigh(coefficients, tolerance=1e-5, verbose=False)`: Get detailed Rayleigh check

### Lorentzian Functions

- `is_lorentzian(coefficients, tolerance=1e-5)`: Check if a polynomial is Lorentzian
- `check_lorentzian(coefficients, tolerance=1e-5, verbose=False)`: Get detailed Lorentzian check

### Utility Functions

- `polynomial_roots(coefficients, tolerance=1e-5)`: Compute polynomial roots
- `evaluate_polynomial(coefficients, x)`: Evaluate polynomial at given value(s)

## Coefficient Format

All functions expect polynomial coefficients in **descending order of powers**:
- `[1, -3, 2]` represents `x^2 - 3x + 2`
- `[1, 0, -1]` represents `x^2 - 1`

## License

MIT License

