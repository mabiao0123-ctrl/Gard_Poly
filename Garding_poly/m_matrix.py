import sympy as sp
import random
import Garding as gd
import numpy as np
def random_nonneg_matrix(n, min_val=0, max_val=100,sparsity=1):
    """
    Generate a random n x n matrix with non-negative integer entries between min_val and max_val. sparsity is the probability of an entry being non-zero.
    """
    return sp.Matrix([[random.randint(min_val, max_val) if random.random()<=sparsity else 0 for _ in range(n)] for _ in range(n)])
    #return sp.Matrix([[random.randint(min_val, max_val) for _ in range(n)] for _ in range(n)])

def is_Z_matrix(A):
    """
    Check if a matrix is a Z-matrix (i.e., all off-diagonal entries are non-positive).
    """
    n = A.shape[0]
    for i in range(n):
        for j in range(n):
            if i != j and A[i, j] > 0:
                return False
    return True

def is_M_matrix(A):
    """
    Check if a matrix is an M-matrix (i.e., a Z-matrix with all eigenvalues having non-negative real parts).
    """
    n = A.shape[0]
    if not is_Z_matrix(A):
        return False
    # Check all eigenvalues are real and non-negative
    eigvals = A.eigenvals()
    for ev in eigvals:
        ev_val = sp.re(ev.evalf())
        if ev_val < 0:
            return False
    return True

def power_iteration_spectral_radius(A, num_simulations: int = 100):
    # Convert sympy Matrix to numpy array of floats
    if isinstance(A, sp.Matrix):
        A = np.array(A).astype(np.float64)
    b_k = np.random.rand(A.shape[1])
    for _ in range(num_simulations):
        b_k1 = np.dot(A, b_k)
        b_k1_norm = np.linalg.norm(b_k1)
        if b_k1_norm == 0:
            # All-zero vector, spectral radius is zero
            return 0.0
        b_k = b_k1 / b_k1_norm
    return np.linalg.norm(np.dot(A, b_k))

def multivariate_characteristic_polynomial(A):
    """
    Returns the multivariate characteristic polynomial det(diag(x1,...,xn) +A)
    for a given sympy Matrix A.
    """
    n = A.shape[0]
    lambdas = sp.symbols('x1:%d' % (n+1))  # x1, x2, ..., xn
    diag_lambda = sp.diag(*lambdas)
    char_poly_multi = (diag_lambda + A).det(method='berkowitz')
    #numer, denom = sp.fraction(char_poly_multi)
    poly = sp.Poly(char_poly_multi, *lambdas)
    return poly

def generate_random_M_matrix(n, min_val=0, max_val=100,sparsity=1):
    A = random_nonneg_matrix(n, min_val, max_val,sparsity=sparsity)
    spectral_radius = power_iteration_spectral_radius(A)
    s=random.randint(int(spectral_radius),2*int(spectral_radius))+1
    M = (s * sp.eye(n) - A)
    return M

def homog_inverse_multi_var_char_poly(A):
    n = A.shape[0]
    lambdas = sp.symbols('x0:%d' % (n+1))  # x1, x2, ..., xn
    diag_lambda = sp.diag(*lambdas[1:])
    char_poly_multi = (lambdas[0]*sp.eye(n)+ diag_lambda*A).det(method='berkowitz')
    #numer, denom = sp.fraction(char_poly_multi)
    poly = sp.Poly(char_poly_multi, *lambdas)
    return poly
def generate_random_M_matrix(n, min_val=0, max_val=100,sparsity=1):
    A = random_nonneg_matrix(n, min_val, max_val,sparsity=sparsity)
    spectral_radius = power_iteration_spectral_radius(A)
    s=random.randint(int(spectral_radius),2*int(spectral_radius))+1
    M = (s * sp.eye(n) - A)
    return M
def homogenized_multivariate_characteristic_polynomial(A):
    """
    Returns the multivariate characteristic polynomial det(diag(x1,...,xn) - A)
    for a given sympy Matrix A.
    """
    n = A.shape[0]
    lambdas = sp.symbols('x0:%d' % (n+1))  # x1, x2, ..., xn
    diag_lambda = sp.diag(*lambdas[1:])
    char_poly_multi = (diag_lambda +lambdas[0]*A).det(method='berkowitz')
    #numer, denom = sp.fraction(char_poly_multi)
    poly = sp.Poly(char_poly_multi, *lambdas)
    return poly
