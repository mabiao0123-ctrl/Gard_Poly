# Represent a multivariate polynomial using monomial basis
# For variables x, y, z and degree up to 3
# P(x, y, z) = 2*x**3 + 5*x*y**2*z + 7*y + 1

import sympy as sp
import itertools
import numpy as np

def is_M_convex(A):
    """
    Check if a finite set A (subset of N^n) is M-convex.
    A should be a list of tuples, each tuple of length n.
    """
    A_set = set(A)
    n = len(A[0])
    for alpha in A:
        for beta in A:
            for i in range(n):
                if alpha[i] > beta[i]:
                    # Find j such that alpha_j < beta_j and alpha - e_i + e_j in A
                    found = False
                    for j in range(n):
                        if alpha[j] < beta[j]:
                            alpha_new = list(alpha)
                            alpha_new[i] -= 1
                            alpha_new[j] += 1
                            if min(alpha_new) < 0:
                                continue
                            if tuple(alpha_new) in A_set:
                                found = True
                                break
                    if not found:
                        return False
    return True

# print("Sympy polynomial:", poly_expr)
# n = 4  # number of variables
# d = 3  # degree
# vars = sp.symbols('x0:%d' % n)  # (x0, x1, x2, x3)

# # Example: monomial basis dictionary for n variables, degree <= d
# poly_dict = {
#     (3, 0, 0, 0): 8,
#     (1, 2, 0, 0): 5,
#     (0, 1, 1, 1): 7,
#     (0, 0, 0, 0): 1,
#     (2, 1, 0, 0): 0,
#     (0, 0, 2, 1): 0,
# }

# # Build sympy polynomial from monomial basis
# poly_expr = sum(
#     coeff * sp.Mul(*[v**e for v, e in zip(vars, exps)])
#     for exps, coeff in poly_dict.items()
# )
# print("Sympy polynomial:", poly_expr)
# def support(poly):
#     """
#     Support of a polynomial: the set of variables with non-zero coefficients.
#     """
#     return {exps for exps, coeff in poly.items() if coeff != 0}
#print("Support of the polynomial:", support(poly_dict))
# #print("Support of the polynomial:", poly.support())

def hessian(expr, variables):
    """
    Compute the Hessian matrix of a sympy expression with respect to the given variables.
    Returns a sympy Matrix of second derivatives.
    """
    n = len(variables)
    H = sp.Matrix(n, n, lambda i, j: sp.diff(expr, variables[i], variables[j]))
    return H

def matrix_eigenvalues(matrix):
    """
    Compute the eigenvalues of a sympy Matrix.
    Returns a dictionary {eigenvalue: algebraic multiplicity}.
    """
    return matrix.eigenvals()

# Example usage:

def count_positive_eigenvalues(matrix):
    """
    Count the number of positive (real) eigenvalues of a sympy Matrix.
    """
    eigenvals = matrix.eigenvals()
    count = 0
    for eig in eigenvals:
        # Try to evaluate the eigenvalue numerically
        try:
            val = sp.N(eig)
            if abs(sp.im(val)) < 1e-10 and val > 0:  # Consider small imaginary part as zero
                count += eigenvals[eig]
        except Exception:
            pass  # If cannot evaluate, skip
    return count


def all_multi_indices(n, k):
    """
    Generate all n-tuples of non-negative integers summing to k.
    """
    for c in itertools.combinations_with_replacement(range(n), k):
        yield c
    

def check_partial_hessians(poly, variables=None):
    """
    Given a degree d polynomial (as sympy expr), check all (d-2)-th order partial derivatives' Hessians.
    For each such derivative, compute its Hessian and check if it has at most 1 positive eigenvalue.
    Return True if all pass, else False.
    """
    poly, expr, poly_dict = poly_pretreat(poly, variables)
    d = sp.total_degree(expr)
    n = len(variables)
    if d < 2:
        return True  # No (d-2)-th derivatives to check
    for idxs in all_multi_indices(n, d-2):
        # Take (d-2)-th partial derivative
        deriv = expr
        for var_idx in idxs:
            
            deriv = sp.diff(deriv, variables[var_idx])
            if deriv == 0:
                    break
        if deriv == 0:
            continue
        # Compute Hessian
        H = hessian(deriv, variables)
        # Count positive eigenvalues
        num_pos = count_positive_eigenvalues(H)
        if num_pos > 1:
            return False
    return True

# Example usage:
# result = check_partial_hessians(poly, vars)
# print("All (d-2)-th partial Hessians have at most 1 positive eigenvalue:", result)
def support(poly, variables=None):
    """
    Support of a polynomial: the set of variables with non-zero coefficients.
    """
    poly, expr, poly_dict = poly_pretreat(poly, variables)
    return {exps for exps, coeff in poly_dict.items() if coeff != 0}

def poly_pretreat(poly, variables=None):
    """
    Preprocess a polynomial to ensure it is in a suitable form for Lorentzian checks.
    Converts sympy expression to a dictionary representation.
    """
    expr=None
    poly_dict = {}
    if isinstance(poly, sp.Expr):
        expr = poly
        poly = sp.Poly(poly, *variables)
        poly_dict = sp.Poly(poly, *variables).as_dict()
    elif isinstance(poly, sp.Poly):
        expr = poly.as_expr()
        poly_dict = poly.as_dict()
        # Already in Poly form
    elif isinstance(poly, dict):
        poly_dict = poly
        if variables is None:
            variables = sp.symbols('x0:%d' % len(next(iter(poly_dict))))  
        poly = sp.Poly.from_dict(poly_dict, variables)
        expr = poly.as_expr()
    return poly, expr, poly_dict

def is_Lorentzian(poly, variables=None):
    """
    Check if a polynomial is Lorentzian:
    1. All coefficients are nonnegative.
    2. The support is M-convex.
    3. All (d-2)-th order partial derivatives have Hessians with at most one positive eigenvalue.
    """
    poly, expr, poly_dict = poly_pretreat(poly, variables)
    #make poly polynomial, expr sp.expr value and poly_dict dict value.
    # Check nonnegative coefficients
    if any(coeff < 0 for coeff in poly_dict.values()):
        return False
    
    # Check M-convexity of support
    supp = list(support(poly_dict))
    if not is_M_convex(supp):
        return False
    
    # Check Hessian condition
    if not check_partial_hessians(expr, variables):
        return False
    
    return True
# x,y,z= sp.symbols('x y z')
# f1=sp.Poly(x**3+2*x**2*y+2*x*y**2, x, y)
# f1.as_expr()
# print("Polynomial:", f1.as_expr())
# print("Is Lorentzian:", is_Lorentzian(f1, (x, y)))

# for index in  all_multi_indices(5, 3):
#     print("Multi-index:", index)
def elementary_symmetric_poly(n, d, variables=None):
    """
    Return the elementary symmetric polynomial of degree d in n variables as a sympy expression.
    """
    if variables is None:
        variables = sp.symbols('x0:%d' % n)
    terms = []
    for idxs in itertools.combinations(range(n), d):
        term = 1
        for i in idxs:
            term *= variables[i]
        terms.append(term)
    return sum(terms)

# Example usage:
# n, d = 10, 8
# vars = sp.symbols('x0:%d' % n)
# esp = elementary_symmetric_poly(n, d, vars)
# print(f"Elementary symmetric polynomial of degree {d} in {n} variables:", esp)
# print("Is Lorentzian:", is_Lorentzian(esp, vars))