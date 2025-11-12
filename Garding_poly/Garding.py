# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import math
import matplotlib.pyplot as plt
# the following compute the sigma_k of lam
def sigma(k,lam):
    n=len(lam)
    if n<k:
        return 0
    if k<0:
        return 0
    if k==0:
        return 1
    if k==1:
        s=0
        for i in range(n):
            s=s+lam[i]
        return s
    s_1=sigma(k-1,lam[:n-1])
    s_2=sigma(k,lam[:n-1])
    return lam[n-1]*s_1+s_2

# compute "sigma_k(lambda:[ij])" by sigma_e(k,lambda,[ij])
def sigma_e(k,lam,excp):
    lam_2=lam.copy()
    for i in excp:
        lam_2[i-1]=0
    return sigma(k,lam_2)

# evaluating Sum_i coeff[n-i] sigma_i(lam)
def fx(coeff,lam):
    s=0
    n=len(coeff)-1
    for i in range(n+1):
        s=s+coeff[i]*sigma(n-i,lam)
    return s

# compute derivative with index in excp
def Dfx(coeff,lam,excp):
    s=0
    j=len(excp)
    l=len(coeff)-1
    for i in range(l+1):
        s=s+coeff[i]*sigma_e(l-i-j,lam,excp)
    return s
        
# compute Hessian matrix for f(lambda)=Sum_i coeff[n-i] sigma_i(lam) and return hessian matrix
def Hess(coeff,la):
    la1=la.copy()
    n=len(la1)
    sw= np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            if i==j:
                sw[i][j]=0
            else:
                sw[i][j]=Dfx(coeff,la1,[i+1,j+1])
    return sw
    
# compute diagonalized f(x) and m-th derivative 
# x=(x_1,...,x_n) the point is each c_k sigma_k(x)=c_k comb(n,k)t^k if x_i=t for all i
# return coefficients
def diagf(coeff,n,m=0):
    "compute diagonalized f(x) and m-th derivative. x=(x_1,...,x_n) the point is each c_k sigma_k(x)=c_k comb(n,k)t^k if x_i=t for all i.  return coefficients"
    c=coeff.copy()#[:-m]
    d=len(c)-1
    for i in range(d+1-m):
        c[i]=math.comb(n,d-i)*math.comb(d-i,m)*c[i]
    return c[:d-m+1]


def diagf1(coeff,m=0):
    "compute m-th derivative of f(t) the difference from diagf is f is really a univariate polynomial"
    c=coeff.copy()#[:-m]
    d=len(c)-1
    for i in range(d+1-m):
        c[i]=math.comb(d-i,m)*c[i]
    return c[:d-m+1]

# dervf(coeff) compute biggest roots of univariate f(x) and derivatives with coeffecients and return a list which consists of s_0,s_1,s_2... where s_i is the biggest root of i-th derivative of f
def dervf(coeff):
    "dervf(coeff) compute biggest roots of univariate f(x) and derivatives with coeffecients and return a list which consists of s_0,s_1,s_2... where s_i is the biggest root of i-th derivative of f"
    s=[]
    d=len(coeff)-1
    for i in range(d):
        coefs=diagf1(coeff,i)
        r=np.roots(coefs)
        real_roots = list(r.real[abs(r.imag)<1e-8]) # where I chose 1-e5 as a threshold
        real_roots.sort(reverse=True)
        if real_roots==[]:
            return "no way!"
        else:
            s.append(real_roots[0].round(5))
    return s


# compute biggest roots of diagonalized f(x) and derivatives with coeffecients 
# return a list which consists of s_0,s_1,s_2... where s_i is the biggest root of i-th derivative of f
# if f has no real roots, return "no way!"
def dervdiagf(coeff,n):
    "compute biggest roots of diagonalized f(x) and derivatives with coeffecients.  return a list which consists of s_0,s_1,s_2... where s_i is the biggest root of i-th derivative of f if f has no real roots, return 'no way!'"
    s=[]
    d=len(coeff)-1
    flag=True
    for i in range(d):
        coefs=diagf(coeff,n,i)
        r=np.roots(coefs)
        real_roots = list(r.real[abs(r.imag)<1e-8]) # where I chose 1-e5 as a threshold
        real_roots.sort(reverse=True)
        if real_roots==[]:
            return "no way!"
        else:
            s.append(real_roots[0].round(5))
    return s

def UpsilonStableDiagnol(coeff,n):
    "Check if a symmetric multivariable polynomial is Upsilon stable"
    r=dervdiagf(coeff,n)
    d=len(coeff)-1
    flag=True
    if r=="no way!":
        return False
    
    for i in range(1,d):
        if r[i]>r[i-1]:
            return False
    return flag

def UpsilonStable(coeff):
    "Check if a univariate polynomial is Upsilon stable"
    if coeff[0]<0:
        return False
    r=dervf(coeff)
    d=len(coeff)-1
    flag=True
    if r=="no way!":
        return False
    
    for i in range(1,d):
        if r[i]>r[i-1]:
            return False
    return flag

# evalf calculate univariate polynomial  f at x 
def evalf(coeff,x):
    d=len(coeff)-1
    s=0
    t=float(x)
    for i in range(d+1):
        s=s+coeff[d-i]*(t**i)
    return s

# constr_f constructs from roots a polynomial f(x) such that the f^(i) has root roots[i]. 
# if roots are sorted in decreasing order then constr_f is the inverse of dervf
def constr_f(roots):
    "construct a Upsilon_stable polynomial using a root sequence"
    "roots:list"
    if roots=="no way!":
        print("cannot construct!")
        return 
    r=roots
    d=len(r)
    a=np.zeros((d+1,d+1))
    a[0][0]=1
    for i in range(1,d+1):
        for j in range(1,i+1):
            a[i][j]=a[i-1][j-1]*i/j
        c=a[i][0:i+1]
        a[i,0]=-evalf(c[::-1],r[d-i])
    return a[d][::-1]

# constr_diagf construct polarized multivariate polynomial from roots of a polynomial f(x) such that the f^(i) has root roots[i]
def constr_diagf(roots,n):
    if roots=="no way!":
        print("cannot construct!")
        return 
    r=roots
    d=len(r)
    a=np.zeros((d+1,d+1))
    a[0][0]=1
    for i in range(1,d+1):
        for j in range(1,i+1):
            a[i][j]=a[i-1][j-1]*i/j
        c=a[i][0:i+1]
        a[i,0]=-evalf(c[::-1],r[d-i])
    s=a[d][::-1]
    for i in range(d):
        s[i]=s[i]/math.comb(n,d-i)
    return s


def specialize_ascend(f,x,a):
    "f arry, n integer, x & a are arrays. f is the coeffcients (in ascending order) of a n variable symmetric multi-affine polynomial. Compute specialize_ascendd univariate polynomial at direction a, specialize_ascend(f,x,a)=coefficient(f(x+at)). the coefficient in an ascending order"
    deg=len(f)-1
    p=np.zeros(deg+1)
    p=p.tolist()
    var=len(x)
    coeff=f
    x_replic=x.copy()
    a_replic=a.copy()
    if deg==0:
        return [coeff[0]]
    if var==0:
        return[coeff[0]]
    if var==1:
        p[0]=coeff[0]+coeff[1]*x[0]
        p[1]=coeff[1]*a[0]
        return p
        
    # if len(x) !=n or len(a) !=n:
    #     return "parameter numbers do not match"
    x_replic.pop(0)
    a_replic.pop(0)
    p1=specialize_ascend(f[1:], x_replic, a_replic)
    p1.append(0)
   # p1.append(0)
    p2=specialize_ascend(f, x_replic, a_replic)
    p[0]=p1[0]*x[0]+p2[0]
    
    for i in range(1,deg+1):
        p[i]=p1[i-1]*a[0]+p1[i]*x[0]+p2[i]
    return p

def specialize(f,x,a):
    "f arry, n integer, x & a are arrays. f is the coeffcients (in descending order) of a n variable symmetric multi-affine polynomial. Compute specialize_ascendd univariate polynomial at direction a, specialize_ascend(f,x,a)=coefficient(f(x+at)). the coefficient in descending order"

    f_ascend=f[::-1]
    p=specialize_ascend(f_ascend,x,a)
    return p[::-1]

def realstable_univarite(f):
    "Show that the univariate polynomial is stable or not"
    r=np.roots(f)
    d=len(f)-1
    real_roots = list(r.real[abs(r.imag)<1e-8]) # where I chose 1-e5 as a threshold
    return len(real_roots)==d
    
# coeff_reverse=[60, -88, 19, 17, -9, 1]
# realstable_univarite(coeff_reverse[::-1])
#z=specialize([1,2,1,2,1,2,1,2,1],[0,0,0,0,0,0,0,0,0,0,0,0],[1,1,1,1,1,1,1,1,1,1,1,1])
#print(z)

def restrict_coeff(c, n, m):
    """
    c coefficients, n variable, degree m.  get the polynomial with x1=0, ([4,18,0,1],4,3)->[1,9,0,1], ([1,9,0,1],3,3)->[3,0,1]; 
    """
    coeff=c[::-1]
    l=n-1
    if l>m:
        l=m
    g=(l+1)*[0]
    for i in range(l+1):
        g[i]=coeff[i]*(n-i)/n
    return g[::-1]

def diagonalize(f, n, m):
    """
    f is the coefficients of a symmetric polynomial in descending order.
    n is the number of variables.
    m is the degree of the polynomial.
    Return the coefficients of the diagonalized polynomial.
    """
    coeff = f[::-1]
    l = n - 1
    if l > m:
        l = m
    g = (l + 1) * [0]
    for i in range(l + 1):
        g[i] = coeff[i] * (n - i) / n
    return g[::-1]

def poly_product(coeffs1, coeffs2):
    """
    Given two lists of coefficients (descending order), return the coefficients of their product (also descending order).
    """
    return np.convolve(coeffs1, coeffs2).tolist()

# Represent a multivariate polynomial using monomial basis
# For variables x, y, z and degree up to 3
# P(x, y, z) = 2*x**3 + 5*x*y**2*z + 7*y + 1

import sympy as sp
import itertools

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
    deg=sp.total_degree(expr)
    H=None
    if deg<2:
        H=sp.Matrix(n,n, lambda i,j:0)
    if deg>2:
        H = sp.Matrix(n, n, lambda i, j: sp.diff(expr, variables[i], variables[j]))
    if deg==2:
        poly=sp.Poly(expr,*variables)
        poly_dict=poly.as_dict()
        H=np.zeros((n,n))
        for i in range(n):
            for j in range(n):
                if i==j:
                    if (2*(i),) in poly_dict:
                        H[i][j]=2*poly_dict[(2*(i),)]
                else:
                    if (1*(i)+1*(j),) in poly_dict:
                        H[i][j]=poly_dict[(1*(i)+1*(j),)]
        H=sp.Matrix(H)
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
    M=np.array(matrix).astype(np.float64)
    eigenvals = np.linalg.eigvals(M)
    count = 0
    for eig in eigenvals:
        # Try to evaluate the eigenvalue numerically
        try:
            val = np.real(eig)
            if abs(np.imag(eig)) < 1e-10 and val > 0:  # Consider small imaginary part as zero
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
    
def second_largest_eigenvalue_iterative(A, num_simulations=30, tol=1e-8):
    """
    Find the second largest real eigenvalue of a real symmetric matrix A using the power iteration and deflation.
    """
    # Convert to numpy array if needed
    if not isinstance(A, np.ndarray):
        A = np.array(A, dtype=np.float64)
    n = A.shape[0]
    # 1. Find largest eigenvalue/vector
    b_k = np.random.rand(n)
    for _ in range(num_simulations):
        b_k1 = np.dot(A, b_k)
        b_k1_norm = np.linalg.norm(b_k1)
        if b_k1_norm == 0:
            break
        b_k = b_k1 / b_k1_norm
    v1 = b_k
    lambda1 = np.dot(v1, np.dot(A, v1))
    # 2. Deflate: project out the top eigenvector
    A_deflated = A - lambda1 * np.outer(v1, v1)
    # 3. Power iteration for second eigenvalue
    b_k = np.random.rand(n)
    b_k -= np.dot(b_k, v1) * v1  # orthogonalize to v1
    for _ in range(num_simulations):
        b_k1 = np.dot(A_deflated, b_k)
        b_k1_norm = np.linalg.norm(b_k1)
        if b_k1_norm == 0:
            break
        b_k = b_k1 / b_k1_norm
        # Re-orthogonalize to v1 to avoid drift
        b_k -= np.dot(b_k, v1) * v1
        b_k /= np.linalg.norm(b_k)
    v2 = b_k
    lambda2 = np.dot(v2, np.dot(A, v2))
    return lambda2

def check_partial_hessians(expr, variables=None):
    """
    Given a degree d polynomial (as sympy expr), check all (d-2)-th order partial derivatives' Hessians.
    For each such derivative, compute its Hessian and check if it has at most 1 positive eigenvalue.
    Return True if all pass, else False.
    """
    if variables is None:
        variables = expr.free_symbols
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
        if deriv == 0 or sp.total_degree(deriv)<=1:
            continue
        # Compute Hessian
        H = hessian(deriv, variables)
        # Count positive eigenvalues
        num_pos = count_positive_eigenvalues(H)
        if num_pos > 1:
            return False
    return True

def check_partial_hessians_eigenvalue_iterative(expr, variables=None):
    """Given a degree d polynomial (as sympy expr), check all (d-2)-th order partial derivatives' Hessians.
    For each such derivative, compute its Hessian and check if it has at most 1 positive eigenvalue using iterative method.
    Return True if all pass, else False."
""" 
    if variables is None:
        variables = expr.free_symbols
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
        # Count positive eigenvalues using iterative method
        second_largest = second_largest_eigenvalue_iterative(H)
        if second_largest > 1e-10:  # Consider numerical tolerance
            return False
    return True


# Example usage:
# result = check_partial_hessians(poly, vars)
# print("All (d-2)-th partial Hessians have at most 1 positive eigenvalue:", result)
def support(poly, variables=None):
    """
    Support of a polynomial: the set of variables with non-zero coefficients.
    """
    poly, _, poly_dict = poly_pretreat(poly, variables)
    return list({exps for exps, coeff in poly_dict.items() if coeff != 0})

def poly_pretreat(poly, variables=None):
    """
    Preprocess a polynomial to ensure it is in a suitable form for Lorentzian checks.
    Converts sympy expression to a dictionary representation.
    """
    expr=None
    poly_dict = {}
    if isinstance(poly, sp.Expr):
        expr = poly
        variables = expr.free_symbols
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

def is_Lorentzian(poly, variables=None,M_convex=False):
    """
    Check if a polynomial is Lorentzian:
    1. All coefficients are nonnegative.
    2. The support is M-convex.
    3. All (d-2)-th order partial derivatives have Hessians with at most one positive eigenvalue.
    """
    expr=None
    poly_dict = {}
    if isinstance(poly, sp.Expr):
        expr = poly
        variables = expr.free_symbols
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
    #make poly polynomial, expr sp.expr value and poly_dict dict value.
    # Check nonnegative coefficients
    if any(coeff < 0 for coeff in poly_dict.values()):
        return False
    
    # Check M-convexity of support
    supp = list(support(poly_dict))
    if not M_convex:
        if not is_M_convex(supp):
            return False
    
    # Check sssHessian condition
    if not check_partial_hessians(expr, variables):
        return False
    
    return True
def is_Lorentzian_numerical(poly, variables=None,M_convex=False):
    """
    Check if a polynomial is Lorentzian:
    1. All coefficients are nonnegative.
    2. The support is M-convex.
    3. All (d-2)-th order partial derivatives have Hessians with at most one positive eigenvalue.
    """
    poly_dict={}
    expr=None   
    if isinstance(poly, sp.Expr):
        expr = poly
        variables = expr.free_symbols
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
    #make poly polynomial, expr sp.expr value and poly_dict dict value.
    # Check nonnegative coefficients
    if any(coeff < 0 for coeff in poly_dict.values()):
        return False
    
    # Check M-convexity of support
    supp = list(support(poly_dict))
    if not M_convex:
        if not is_M_convex(supp):
            return False
    
    # Check sssHessian condition
    if not check_partial_hessians_eigenvalue_iterative(expr, variables):
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

def homogenize_poly(poly, n, d, variables=None):
    """
    Homogenize a polynomial to degree d in n variables by adding a new variable x0.
    """
    poly_dict = {}
    if variables is None:
        variables = sp.symbols('x0:%d' % (n+1))
    x0 = variables[0]
    new_vars = variables[1:]
    new_poly_dict = {}
    if isinstance(poly, sp.Expr):
        poly = sp.Poly(poly, *new_vars)
        poly_dict = poly.as_dict()
    elif isinstance(poly, sp.Poly):
        poly_dict = poly.as_dict()
        # Already in Poly form
    elif isinstance(poly, dict):
        poly_dict = poly
    
    for exps, coeff in poly_dict.items():
        total_deg = sum(exps)
        if total_deg > d:
            raise ValueError("Polynomial degree exceeds target degree for homogenization.")
        new_exps = (d - total_deg,) + exps  # Add exponent for x0
        new_poly_dict[new_exps] = coeff
    return sp.Poly.from_dict(new_poly_dict, (x0,) + new_vars)
# Example usage:
# n, d = 10, 8
# vars = sp.symbols('x0:%d' % n)
# esp = elementary_symmetric_poly(n, d, vars)
# print(f"Elementary symmetric polynomial of degree {d} in {n} variables:", esp)
# print("Is Lorentzian:", is_Lorentzian(esp, vars))

def transform_from_symbol_poly(poly, kappa, variables):
    """"
    Given a symbol polynomial poly in variables, construct a transformation 
    """
    poly_dict={}
    import math
    if isinstance(poly, sp.Expr):
        expr = poly
        poly_dict=sp.Poly(expr,variables).as_dict()
    elif isinstance(poly, sp.Poly):
        poly_dict=poly.as_dict()
    elif isinstance(poly, dict):
        poly_dict=poly
    n=len(kappa)
    m=len(variables)-n
    mapping = {}
    for key in poly_dict:
        var1 = tuple(kappa[i] - key[i] for i in range(n))
        var2 = key[n:]
        adjust_coeff = 1
        for i in range(n):
            if  var1[i]<0:
                print("kappa is smaller than exponent. cannot construct transform")
                return None
            adjust_coeff *= math.factorial(kappa[i]) / math.factorial(var1[i])/ math.factorial(key[i])
        
        coeff = poly_dict[key] / adjust_coeff
        if var1 not in mapping:
            mapping[var1] = {}
        mapping[var1][var2] = coeff
    return mapping

def evaluate_transform_poly(mapping,poly,variables=None):
    """
    Given a transform T and a polynomial f, evaluate T(f).
    """
    poly_dict={}
    if isinstance(poly, sp.Expr):
        expr = poly
        if variables is None:
            variables = list(expr.free_symbols)
        poly_dict=sp.Poly(expr,variables).as_dict()
    elif isinstance(poly, sp.Poly):
        poly_dict=poly.as_dict()
        variables=list(poly.gens)
    elif isinstance(poly, dict):
        poly_dict=poly
        if variables is None:
            variables=sp.symbols('x0:%d' % (len(list(poly_dict.keys())[0])))
        else:
            variables=list(variables) if not isinstance(variables, list) else variables

    poly_dict_new={}
    m=0
    for key in poly_dict:
        if not key in mapping:
            continue
        target=mapping[key]
        
        if variables is None or len(key) != len(variables):
            return print("variable numbers are different. cannot evaluate")
        for support in target:
            m=len(support)
            if support in poly_dict_new:
                poly_dict_new[support]+=poly_dict[key]*target[support]
            else:
                poly_dict_new[support]=target[support]*poly_dict[key]
    vars=sp.symbols('x0:%d' % m)
    return sp.Poly.from_dict(poly_dict_new,vars)

def plot_poly(poly, plot_range=10, points=300):
    """Plot the polynomial and its positive/negative regions; accepts a sympy Expr or Poly in two variables."""
    # Ensure we have a sympy Poly object
    if isinstance(poly, sp.Expr):
        poly = sp.Poly(poly)
    else:
        try:
            poly = sp.Poly(poly)
        except Exception:
            raise ValueError("plot_poly: input must be a sympy Expr or Poly convertible to Poly")

    if not isinstance(poly, sp.Poly):
        raise ValueError("plot_poly: input must be a sympy Expr or Poly")

    gens = poly.gens
    if len(gens) != 2:
        raise ValueError("plot_poly requires a polynomial in exactly two variables")
    x, y = gens

    f_np = sp.lambdify((x, y), poly.as_expr(), 'numpy')
    # Create a grid
    X, Y = np.meshgrid(
        np.linspace(-plot_range, plot_range, points),
        np.linspace(-plot_range, plot_range, points),
    )
    Z = f_np(X, Y)

    plt.figure(figsize=(6, 6))
    # Plot the zero set (contour where f=0)
    plt.contour(X, Y, Z, levels=[0], colors='k', linewidths=2)

    # Compute bounds for filled contours, guarding against degenerate cases
    zmin = float(np.nanmin(Z))
    zmax = float(np.nanmax(Z))
    if zmin == 0:
        zmin -= 1e-8
    if zmax == 0:
        zmax += 1e-8

    neg_levels = [zmin, 0]
    pos_levels = [0, zmax]

    # Fill f<0 region
    plt.contourf(X, Y, Z, levels=neg_levels, colors=['#f38ba8'], alpha=0.5)

    # Fill f>0 region
    plt.contourf(X, Y, Z, levels=pos_levels, colors=['#a6e3a1'], alpha=0.5)

    plt.xlabel(str(x))
    plt.ylabel(str(y))
    plt.title('Regions where f(x, y) > 0 (green) and f(x, y) < 0 (red)')
    plt.grid(True)
    plt.show()
    return None

def homogenization(poly,deg=None,vars=None,Z=None):
    "Accept polynomials from sympy poly class and return an homogeneous polynomial with additional variable Z"
    "poly: sympy Poly class"
    "deg: target degree of the homogeneous polynomial; if None, use poly's total degree"
    "vars: variables to use for homogenization; if None, use poly's gens"
    "Z: symbol for the homogenizing variable; if None, use 'Z'"
    
    if deg is None:
        deg = poly.total_degree()
    # normalize vars to a tuple of symbols; allow single Symbol or iterable
    if vars is None:
        vars = poly.gens()
    else:
        if isinstance(vars, sp.Symbol):
            vars = (vars,)
        else:
            try:
                len(vars)
            except TypeError:
                vars = (vars,)
        # if the provided vars length doesn't match the polynomial's gens, fall back to poly.gens()
    # ensure vars is a tuple so concatenation with (Z,) works and indexing is stable
    vars = tuple(vars)
    if Z is None:
        Z=sp.symbols('Z')
    homog_f = 0
    for monom, coeff in poly.terms():
        monom_degree = sum(monom)
        Z_degree = deg - monom_degree
        homog_monom = coeff
        for i in range(len(vars)):
            var = vars[i]
            var_degree = monom[i]
            homog_monom *= var**var_degree
        homog_monom *= Z**Z_degree
        homog_f += homog_monom
    return sp.Poly(homog_f, vars + (Z,))

def check_non_neg_num_rand(f,x,R=1,num_samples=10000):
    "Check if a function is positive for all variables in [0,R]"
    "f: sympy polynomial (Poly or Expr)"
    "x: list of sympy symbols representing variables"
    "R: upper bound of the box [0,R]^n"
    "num_samples: number of random samples to test"
    if isinstance(f, sp.Poly):
        f = f.as_expr()
    
    from sympy import lambdify
    f_num = lambdify(x, f, 'numpy')
    n_vars = len(x)
    all_nonneg = True
    for _ in range(num_samples):
        point = np.random.uniform(0, R, n_vars)
        val = f_num(*point)
        if val < 0:
            return False,point,val
    return True,0,0

def check_non_neg_opt(f, vars, R=1):
    from scipy.optimize import minimize
    from sympy import lambdify
    # convert sympy.Poly to sympy.Expr (call the method, not reference it)
    if isinstance(f, sp.Poly):
        f = f.as_expr()
    # Create a numeric function from sympy expression.
    # `vars` should be an iterable of sympy symbols, e.g. [x, y, z]
    f_num = lambdify(vars, f, 'numpy')

    n_vars = len(vars)

    # bounds enforce 0 <= xi <= R
    bounds = [(0, R) for _ in range(n_vars)]

    # constraints need callables that accept the optimizer vector `v`
    constraints = []
    for i in range(n_vars):
        # R - v[i] >= 0  -> v[i] <= R
        constraints.append({'type': 'ineq', 'fun': (lambda i: (lambda v: R - v[i]))(i)})
        # v[i] >= 0 is already enforced by bounds, so no extra constraint needed

    # Wrapper so the objective accepts a single vector argument
    def obj(v):
        # ensure f_num is called with separate args
        try:
            val = f_num(*v)
        except TypeError:
            # if f_num expects a single sequence-like
            val = f_num(v)
        return float(val)

    # Try several random starting points to avoid local minima
    min_val = np.inf
    pts=None
    for _ in range(20):
        x0 = np.random.uniform(0, R, n_vars)
        res = minimize(obj, x0, bounds=bounds, constraints=constraints, method='SLSQP')
        if res.success and res.fun < min_val:
            min_val = res.fun
            pts=res.x

    if min_val >= 0:
        return True, min_val, pts
    else:
        return False, min_val, pts
    
def check_Rayleigh(f, vars, Method="opt", trials=10000,option="show_arg"):
    
    "check if a polynomial f satisfies Rayleigh's condition using numerical or optimization methods"
    "f: sympy polynomial (Poly or Expr)"
    "vars: list of sympy symbols representing variables"
    "Method: 'num' for numerical sampling, 'opt' for optimization"
    "trials: number of samples for 'num' method only required for 'num' method"
    "option: 'show_arg' to print counterexample arguments if found, else no print"
     
    n = len(vars)
    # Ensure f is a sympy.Poly over the provided variables.
    # Always use sp.Poly (avoid unqualified Poly to prevent NameError).
    if not isinstance(f, sp.Poly):
        try:
            # pass the generator list as a sequence (safer than star-unpacking)
            f = sp.Poly(f, vars)
        except Exception:
            # fallback if f is a tuple like (poly, gens)
            if isinstance(f, tuple) and len(f) > 0:
                f = sp.Poly(f[0], vars)
            else:
                if isinstance(f,sp.Expr):
                    f=sp.Poly(f,vars)
                else:
                    raise

    def is_homogeneous(poly, vars):
        d = poly.total_degree()
        f_dict = poly.as_dict()
        for expr in f_dict:
            if sum(expr) != d:
                return False
        return True
    if not is_homogeneous(f, vars):
        f=homogenization(f,deg=f.total_degree(),vars=vars,Z=None)
        vars=list(f.gens)
    for i in range(n):
        for j in range(i, n):
            # Work with polynomial expressions (as sympy.Expr) for differentiation,
            # then convert back to Poly if needed.
            fi_expr = f.as_expr().diff(vars[i])
            fj_expr = f.as_expr().diff(vars[j])
            fij_expr = fi_expr.diff(vars[j])

            # build Wronskian as sympy expression
            Wronskian = f.as_expr() * fij_expr - fi_expr * fj_expr

            if Method == "num":
                is_nonneg, point, val = check_non_neg_num_rand(-1 * Wronskian, vars, R=1, num_samples=trials)
                if not is_nonneg and option=="show_arg":
                    print(is_nonneg, point, val)
                    print(Wronskian,"i=",i,"j=",j)
                return is_nonneg
            elif Method == "opt":
                is_nonneg, point, val = check_non_neg_opt(-1 * Wronskian, vars, R=1)
                if not is_nonneg and option=="show_arg":
                    print(is_nonneg, point, val)
                    print(Wronskian,"i=",i,"j=",j)
                return is_nonneg
    return True