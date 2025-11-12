#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 25 14:09:50 2025

@author: biaoma
"""
import random
import sigma
import numpy as np
#from sympy import *
#x,y,z,w,t = symbols("x,y,z,w,t")
#p= x**2 +1+y
#p.all_coeffs()
# r=[4,3,2,1,0]
# c=sigma.constr_f(r)
# sigma.dervf(c)

# c=[-7.42691, -210.933, -164.909, 0, 1]
# c=c[::-1]
# print(c)
# sigma.dervf(c)
#p=p.subs([(x,3+5*t),(y,2+5*t)])#.all_coeffs()
#Poly(p,t).all_coeffs()
# random.seed()
# flag="Yes"
# for k in range(1000):
    
#      F=x*y*z*w -(random.uniform(0,10) *x*y+random.uniform(0,100) *y *z +random.uniform(0,8) *z* w +random.uniform(0,21) *w* x + random.uniform(0,51)* y* w+random.uniform(0,12)* x* z) - (random.uniform(0,13)* x +random.uniform(0,166)* y  +random.uniform(0,82)* z  +random.uniform(0,21)* w  )-random.uniform(0,10)
#      F=F.subs([(x,random.uniform(-10,10)+random.uniform(0,10)*t),(y,random.uniform(-10,10)+random.uniform(0,10)*t),(z,random.uniform(-10,10)+random.uniform(0,10)*t),(w,random.uniform(-10,10)+random.uniform(0,10)*t)])
#      #print(Poly(F,t))
#      coeff=Poly(F,t).all_coeffs()
#      if sigma.UpsilonStable(coeff)==False:
#          flag=False
# print(flag)



# def proj(x,excp):
#     "project x to the hyperplane  x_i=0 for i in excp"
#     p=x
#     p.pop(excp-1)
#     return p


# f=sigma.constr_diagf([4,3,2,1], 4)
# #print(np.zeros(3))
# a=sigma.specialize(f,[1,3,1,-200],[1,300,0.1,0.5])
# b1=sigma.dervdiagf(f, 4)
# b=sigma.UpsilonStableDiagnol(f,4)
# c=sigma.UpsilonStable(a)

random.seed()
n=6 #degree of polynomial
m=n+1 # variable number x_1,x_2,...,x_m
count_realstable=0 # couning real stable polynomials
flag=True
for l in range(100):
    
    # generate roots of derivatives in descending order
    roots=[0] # roots start with 0
    for i in range(1,n):
        new=roots[i-1]-random.uniform(0.2,15)**3
        roots.append(new)
    
    # construct polynomial by roots of derivatives
    f=sigma.constr_diagf(roots, m)
    
    # exclude real stable case which is trivial
    if  sigma.realstable_univarite(f):
        count_realstable+=1
    else:
        stable=True
        for s in range(5):
            x=[]
            for i in range(m):
                x.append(random.uniform(-100000,100000))
            for k in range(50):
                a=[]
                for i in range(m):
                    a.append(random.uniform(0,10))
                f_spec=sigma.specialize(f, x, a)
                r=sigma.dervf(f_spec)
                stable1= sigma.UpsilonStable(f_spec)
                if stable1==False:
                    print("x=",x)
                    print("a=",a)
                    print("f_spec=",f_spec)
                    print("f_spec roots=", r)
                    stable=stable and stable1        
                #print("f_spec is Upsilon stable")
        flag=flag and stable
if flag:
    print("all polynomials are Upsilon stable")
print("numbers of real stable polynomials", count_realstable)
