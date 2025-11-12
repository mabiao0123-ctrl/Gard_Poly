import math
import random
import sigma
import numpy as np
def taylor_coeffs_at_t0(coeffs, t0):
    """
    Given coeffs [a0, a1, ..., an] for f(t) = a0 + a1*t + ... + an*t^n,
    return [c0, c1, ..., cn] for f(t) = c0 + c1*(t-t0) + ... + cn*(t-t0)^n
    """
    n = len(coeffs) - 1
    new_coeffs = []
    for j in range(n+1):
        cj = 0
        for k in range(j, n+1):
            cj += coeffs[k] * math.comb(k, j) * (t0**(k-j))
        new_coeffs.append(cj)
    return new_coeffs

p=taylor_coeffs_at_t0([1,4,6,4,1], -1)
print(p)  # Output: [1, 1, 0] for f(t) = 1 + 2*t + t^2 at t0 = 1
random.seed()
n=4 #degree of polynomial
m=n # variable number x_1,x_2,...,x_m
count_realstable=0 # couning real stable polynomials
flag=True
# generate roots of derivatives in descending order
roots=[-1,-2,-2] # roots start with 0

# construct polynomial by roots of derivatives
print("roots=", roots)
# roots=[-1,-2,-3,-4]
f=sigma.constr
print("f=", f)
f_inverse=f[::-1]
print("f_inverse is Upsilon stable:", sigma.UpsilonStable([5,6,12,1]))
print("roots of f_inverse", sigma.dervf(f_inverse))