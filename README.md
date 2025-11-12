# Gard_Poly

A Python package for testing and analyzing special classes of polynomials including real stable, Lorentzian, Rayleigh, and Gårding polynomials.

## Table of Contents
- [Features](#features)
- [Installation](#installation)
- [Dependencies](#dependencies)
- [Quick Start](#quick-start)
- [Main Functions](#main-functions)
- [Examples](#examples)
- [Theory Background](#theory-background)
- [Contributing](#contributing)

## Features

- **Polynomial Classification**: Test polynomials for various stability properties
  - Real stable polynomials
  - Lorentzian polynomials
  - Rayleigh polynomials
  - $\Upsilon$-stable (Upsilon-stable) polynomials (Univariate Garding)
  - M-convex support checking
  
- **Matrix Analysis**: 
  - M-matrix and Z-matrix testing
  - Multivariate characteristic polynomials
  - Spectral radius computation
  - Random M-matrix generation
  
- **Visualization**: Plot polynomial regions and zero sets

- **Analysis Tools**:
  - Compute derivatives and partial derivatives
  - Hessian matrix computation
  - Root finding for univariate and multivariate polynomials
  - Polynomial transformations and homogenization

## Installation

Clone this repository:

```bash
git clone https://github.com/mabiao0123-ctrl/Gard_Poly.git
cd Gard_Poly
```

## Dependencies

Required packages:
- **SymPy**: Symbolic mathematics
- **NumPy**: Numerical computing
- **SciPy**: Scientific computing (optional, for advanced features)
- **Matplotlib**: Plotting and visualization

Install dependencies:

```bash
pip install sympy numpy scipy matplotlib
```

## Quick Start

```python
import sympy as sp
from Garding_poly.Garding import (
    is_Lorentzian, check_Rayleigh, UpsilonStable, dervf,
    elementary_symmetric_poly, plot_poly
)
from Garding_poly.m_matrix import is_M_matrix, multivariate_characteristic_polynomial

# Define variables
x, y, z = sp.symbols('x y z')

# Create a polynomial
poly = x**3 + 2*x**2*y + 2*x*y**2 + y**3

# Check if it's Lorentzian
result = is_Lorentzian(poly, (x, y))
print(f"Is Lorentzian: {result}")

# Check Υ-stability for univariate polynomial
coeff = [1, -3, 3, -1]  # (x-1)³
is_ustable = UpsilonStable(coeff)
roots = dervf(coeff)
print(f"Polynomial is Υ-stable: {is_ustable}")
print(f"Root sequence: {roots}")

# Check Rayleigh property
rayleigh_poly = x**2 + y**2 + z**2
is_rayleigh = check_Rayleigh(rayleigh_poly, [x, y, z])
print(f"Satisfies Rayleigh property: {is_rayleigh}")

# M-matrix analysis
A = sp.Matrix([[3, -1], [-1, 2]])
print(f"Is M-matrix: {is_M_matrix(A)}")
char_poly = multivariate_characteristic_polynomial(A)
print(f"Characteristic polynomial: {char_poly}")

# Generate elementary symmetric polynomial
esp = elementary_symmetric_poly(n=4, d=2)
print(f"Elementary symmetric polynomial: {esp}")

# Visualize a 2-variable polynomial
plot_poly(poly, plot_range=5)
```

## Main Functions

### Polynomial Testing

#### `is_Lorentzian(poly, variables=None, M_convex=False)`
Check if a polynomial satisfies the Lorentzian property:
1. All coefficients are non-negative
2. Support is M-convex
3. All (d-2)-th order partial derivatives have Hessians with at most one positive eigenvalue

**Parameters:**
- `poly`: SymPy expression, Poly object, or dictionary
- `variables`: List of variables (optional, auto-detected if not provided)

**Returns:** Boolean

#### `UpsilonStable(coeff)` 
Check if a univariate polynomial is Υ-stable (Upsilon-stable). A polynomial is Υ-stable if the largest real roots of its successive derivatives form a decreasing sequence: r₀ ≥ r₁ ≥ r₂ ≥ ... ≥ r_{d-1}.

**Parameters:**
- `coeff`: List of coefficients in descending order [aₙ, aₙ₋₁, ..., a₁, a₀]

**Returns:** Boolean

**Example:** For polynomial x³ - 3x² + 3x - 1, use `coeff = [1, -3, 3, -1]`

#### `UpsilonStableDiagnol(coeff, n)`
Check if a symmetric multivariate polynomial (represented by its diagonal restriction) is Υ-stable.

**Parameters:**
- `coeff`: List of coefficients in descending order
- `n`: Number of variables

**Returns:** Boolean

#### `dervf(coeff)`
Find the largest real root of a univariate polynomial and each of its derivatives. Returns the root sequence [r₀, r₁, ..., r_{d-1}] where rᵢ is the largest root of f⁽ⁱ⁾(x).

**Parameters:**
- `coeff`: List of coefficients in descending order

**Returns:** List of floats (root sequence), or "no way!" if no real roots exist

#### `dervdiagf(coeff, n)`
Find root sequence for the diagonal restriction of a symmetric multivariate polynomial.

**Parameters:**
- `coeff`: List of coefficients in descending order
- `n`: Number of variables

**Returns:** List of floats (root sequence), or "no way!" if no real roots exist

#### `constr_f(roots)`
Construct a univariate Υ-stable polynomial from a decreasing root sequence. Inverse operation of `dervf`.

**Parameters:**
- `roots`: List of roots in decreasing order [r₀, r₁, ..., r_{d-1}]

**Returns:** Array of polynomial coefficients in descending order

#### `constr_diagf(roots, n)`
Construct a symmetric multivariate Υ-stable polynomial from a decreasing root sequence.

**Parameters:**
- `roots`: List of roots in decreasing order
- `n`: Number of variables

**Returns:** Array of polynomial coefficients in descending order

#### `check_Rayleigh(f, vars, Method="opt", trials=10000, option="show_arg")`
Check if a homogeneous polynomial satisfies the Rayleigh property using numerical or optimization methods. The Rayleigh condition requires that for all i, j: `f · ∂²f/∂xᵢ∂xⱼ - ∂f/∂xᵢ · ∂f/∂xⱼ <= 0` on the unit sphere.

**Parameters:**
- `f`: SymPy expression or Poly object (will be homogenized if not already homogeneous)
- `vars`: List of SymPy symbols representing variables
- `Method`: `'opt'` for optimization-based checking, `'num'` for random sampling (default: `'opt'`)
- `trials`: Number of random samples (only used when Method='num', default: 10000)
- `option`: `'show_arg'` to print counterexample if found, else no print (default: `'show_arg'`)

**Returns:** Boolean (True if satisfies Rayleigh property, False otherwise)

### Symmetric Polynomial Functions

#### `sigma(k, lam)`
Compute the k-th elementary symmetric polynomial evaluated at λ.

**Parameters:**
- `k`: Integer, degree of symmetric polynomial
- `lam`: List/array of values

**Returns:** Numeric value

#### `elementary_symmetric_poly(n, d, variables=None)`
Generate the elementary symmetric polynomial of degree d in n variables.

**Returns:** SymPy expression

### Hessian Analysis

#### `hessian(expr, variables)`
Compute the Hessian matrix of a SymPy expression.

**Returns:** SymPy Matrix

#### `check_partial_hessians(expr, variables=None)`
Check if all (d-2)-th order partial derivatives have Hessians with at most 1 positive eigenvalue.

**Returns:** Boolean

### Polynomial Transformations

#### `homogenize_poly(poly, n, d, variables=None)`
Homogenize a polynomial to degree d by adding a new variable x₀.

#### `specialize(f, x, a)`
Compute the univariate polynomial f(x+at) along direction a.

**Parameters:**
- `f`: Coefficient array (descending order)
- `x`: Point coordinates
- `a`: Direction vector

**Returns:** List of coefficients

### Visualization

#### `plot_poly(poly, plot_range=10, points=300)`
Plot a 2-variable polynomial showing positive/negative regions.

**Parameters:**
- `poly`: SymPy expression or Poly object
- `plot_range`: Plotting range (default: 10)
- `points`: Grid resolution (default: 300)

### M-Matrix 

#### `is_Z_matrix(A)`
Check if a matrix is a Z-matrix (all off-diagonal entries are non-positive).

**Parameters:**
- `A`: SymPy Matrix

**Returns:** Boolean

#### `is_M_matrix(A)`
Check if a matrix is an M-matrix (a Z-matrix with all eigenvalues having non-negative real parts).

**Parameters:**
- `A`: SymPy Matrix

**Returns:** Boolean

#### `generate_random_M_matrix(n, min_val=0, max_val=100, sparsity=1)`
Generate a random n×n M-matrix.

**Parameters:**
- `n`: Matrix dimension
- `min_val`: Minimum value for matrix entries (default: 0)
- `max_val`: Maximum value for matrix entries (default: 100)
- `sparsity`: Probability of an entry being non-zero, between 0 and 1 (default: 1)

**Returns:** SymPy Matrix

#### `multivariate_characteristic_polynomial(A)`
Compute the multivariate characteristic polynomial det(diag(x₁,...,xₙ) + A) for a given matrix A.

**Parameters:**
- `A`: SymPy Matrix of size n×n

**Returns:** SymPy Poly in variables x₁, x₂, ..., xₙ

#### `homogenized_multivariate_characteristic_polynomial(A)`
Compute the homogenized multivariate characteristic polynomial det(diag(x₁,...,xₙ) + x₀·A).

**Parameters:**
- `A`: SymPy Matrix of size n×n

**Returns:** SymPy Poly in variables x₀, x₁, ..., xₙ

#### `power_iteration_spectral_radius(A, num_simulations=100)`
Compute the spectral radius (largest eigenvalue magnitude) of a matrix using power iteration.

**Parameters:**
- `A`: SymPy Matrix or NumPy array
- `num_simulations`: Number of iterations (default: 100)

**Returns:** Float (spectral radius)

### Root sequence

#### `dervf(coeff)`
Find the largest real roots of a univariate polynomial and its derivatives.

**Returns:** List of roots [r₀, r₁, r₂, ...] where rᵢ is the largest root of f⁽ⁱ⁾

#### `constr_f(roots)`
Construct a Υ-stable polynomial from a decreasing sequence of roots.

**Returns:** Array of polynomial coefficients

## Examples

### Example 1: Testing Lorentzian Polynomials

```python
import sympy as sp
from Garding_poly.Garding import is_Lorentzian

# Define variables
x, y = sp.symbols('x y')

# Test polynomial
poly = x**3 + 3*x**2*y + 3*x*y**2 + y**3

# Check Lorentzian property
is_lor = is_Lorentzian(poly, (x, y))
print(f"Polynomial is Lorentzian: {is_lor}")
```

### Example 2: Elementary Symmetric Polynomials

```python
from Garding_poly.Garding import elementary_symmetric_poly

# Generate e_3(x₀, x₁, x₂, x₃) = x₀x₁x₂ + x₀x₁x₃ + x₀x₂x₃ + x₁x₂x₃
esp = elementary_symmetric_poly(n=4, d=3)
print(esp)
```

### Example 3: Checking Υ-Stability

```python
from Garding_poly.Garding import UpsilonStable, dervf, constr_f

# Test 1: Check if a polynomial is Υ-stable
# Coefficients in descending order: [a_n, a_{n-1}, ..., a_1, a_0]
coeff1 = [1, -3, 3, -1]  # This is (x-1)³
is_stable1 = UpsilonStable(coeff1)
print(f"Polynomial (x-1)³ is Υ-stable: {is_stable1}")

# Test 2: Another example
coeff2 = [1, 0, -1]  # This is x² - 1 = (x-1)(x+1)
is_stable2 = UpsilonStable(coeff2)
print(f"Polynomial x² - 1 is Υ-stable: {is_stable2}")

# Test 3: Find root sequence
coeff3 = [1, -6, 11, -6]  # This is (x-1)(x-2)(x-3)
roots = dervf(coeff3)
print(f"Root sequence of derivatives: {roots}")
print(f"Is decreasing: {all(roots[i] >= roots[i+1] for i in range(len(roots)-1))}")

# Test 4: Construct a Υ-stable polynomial from root sequence
root_sequence = [3.0, 2.5, 2.0]  # Decreasing sequence
poly_coeff = constr_f(root_sequence)
print(f"Constructed Υ-stable polynomial coefficients: {poly_coeff}")
print(f"Verification - is Υ-stable: {UpsilonStable(poly_coeff)}")
```

### Example 4: Hessian Analysis

```python
import sympy as sp
from Garding_poly.Garding import hessian, check_partial_hessians

x, y, z = sp.symbols('x y z')
poly = x**4 + y**4 + z**4 + 2*x**2*y**2 + 2*y**2*z**2 + 2*x**2*z**2

# Compute Hessian
H = hessian(poly, [x, y, z])
print("Hessian matrix:")
print(H)

# Check partial Hessian condition
result = check_partial_hessians(poly, [x, y, z])
print(f"Passes partial Hessian test: {result}")
```

### Example 5: Checking Rayleigh Property

```python
import sympy as sp
from Garding_poly.Garding import check_Rayleigh

# Define variables
x, y, z = sp.symbols('x y z')

# Test a homogeneous polynomial
poly = x**2 +x*y+ y**2 +x*z+y*z+ z**2

# Check Rayleigh property using optimization method
is_rayleigh = check_Rayleigh(poly, [x, y, z], Method="opt")
print(f"Polynomial satisfies Rayleigh property: {is_rayleigh}")

# Alternatively, use numerical sampling with more trials
is_rayleigh_num = check_Rayleigh(poly, [x, y, z], Method="num", trials=50000)
print(f"Polynomial satisfies Rayleigh property (numerical): {is_rayleigh_num}")
```

### Example 6: M-Matrix Operations

```python
import sympy as sp
from Garding_poly.m_matrix import (
    is_M_matrix, 
    is_Z_matrix, 
    generate_random_M_matrix,
    multivariate_characteristic_polynomial
)

# Create a matrix
A = sp.Matrix([
    [3, -1, 0],
    [-1, 3, -1],
    [0, -1, 2]
])

# Check if it's a Z-matrix
print(f"Is Z-matrix: {is_Z_matrix(A)}")

# Check if it's an M-matrix
print(f"Is M-matrix: {is_M_matrix(A)}")

# Generate a random M-matrix
M = generate_random_M_matrix(n=3, min_val=0, max_val=10, sparsity=0.7)
print(f"Random M-matrix:\n{M}")

# Compute multivariate characteristic polynomial
char_poly = multivariate_characteristic_polynomial(A)
print(f"Characteristic polynomial: {char_poly}")
```

## Theory Background

### Lorentzian Polynomials
A polynomial is **Lorentzian** if:
1. It has non-negative coefficients
2. Its support is M-convex
3. All (d-2)-th order partial derivatives have Hessians with at most one positive eigenvalue

Lorentzian polynomials generalize real stable polynomials and appear in combinatorics, optimization, and tropical geometry.

### $\Upsilon$-Stable Polynomials (Univariate Gårding Polynomials)

A univariate polynomial f(x) of degree d is **Υ-stable** (Upsilon-stable) if the largest real roots of its successive derivatives form a decreasing sequence:

$$r_0 \geq r_1 \geq r_2 \geq \cdots \geq r_{d-1}$$

where $r_i$ is the largest real root of $f^{(i)}(x)$ (the i-th derivative).

#### Key Properties:
- **Root Sequence**: The sequence [r₀, r₁, ..., r_{d-1}] completely characterizes the polynomial's stability
- **Constructive**: Given any decreasing sequence of real numbers, one can construct a unique (up to scaling) Υ-stable polynomial with that root sequence
- **Univariate Gårding**: These are univariate analogues of Gårding hyperbolic polynomials


### M-Convex Sets
A set A ⊆ ℕⁿ is **M-convex** if for all α, β ∈ A and all i with αᵢ > βᵢ, there exists j with αⱼ < βⱼ such that α - eᵢ + eⱼ ∈ A.

M-convexity is a discrete analog of convexity important in discrete optimization.

### Rayleigh Polynomials
A homogeneous polynomial f(x₁, ..., xₙ) satisfies the **Rayleigh property** if for all i, j and all points on the unit sphere:

$$f \cdot \frac{\partial^2 f}{\partial x_i \partial x_j} - \frac{\partial f}{\partial x_i} \cdot \frac{\partial f}{\partial x_j} \leq 0$$

This is equivalent to the Wronskian condition being non-positive. Rayleigh polynomials are related to log-concave functions and appear in the study of eigenvalue inequalities.

### M-Matrices and Z-Matrices

A matrix A is a **Z-matrix** if all its off-diagonal entries are non-positive:
$$a_{ij} \leq 0 \text{ for all } i \neq j$$

A **M-matrix** is a Z-matrix whose eigenvalues all have non-negative real parts. Equivalently, an M-matrix can be written as A = sI - B where s > ρ(B) (ρ is the spectral radius) and B has non-negative entries.

#### Properties of M-matrices:
- All principal minors are positive
- The inverse (if it exists) has non-negative entries
- M-matrices arise naturally in numerical analysis, economics, and optimization

#### Multivariate Characteristic Polynomial:
For an n×n matrix A, the multivariate characteristic polynomial is:
$$p(x_1, \ldots, x_n) = \det(\text{diag}(x_1, \ldots, x_n) + A)$$

This generalizes the classical characteristic polynomial det(λI - A) and provides connections between matrix theory and multivariate polynomial analysis. When A is an M-matrix, this polynomial has special stability properties.

