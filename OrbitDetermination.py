# -*- coding: utf-8 -*-
"""
Created on Sat Dec  2 00:03:16 2017

@author: Adam
"""

# Program Details

"""
This program will take an input of two 3 dimensional vectors defining two positions observed from a spacecraft orbiting a massive object and a time of flight
measurement between them, and find the orbit of said object as a function of the radial and tangential positions and velocities. The algorithm used in this 
program is outline in chapter 5 of 'Fundamentals of Astrodynamics' by Bate, Mueller, and White.
"""

# Libraries

from math import factorial

# Functions

def del_nu()

def C(z,epsilon):
    
    k=1
    C_k0 = 1/2
    C_k1 = C_k0 + ((-z)**k)/factorial(2*k+2)
    
    while abs(C_k1 - C_k0) > epsilon:
        k += 1
        C_k0 = C_k1
        C_k1 += ((-z)**k)/factorial(2*k+2)
    
    return C_k1

def S(z,epsilon):
    
    k=1
    S_k0 = 1/factorial(3)
    S_k1 = S_k0 +((-z)**k)/factorial(2*k+3)
    
    while abs(S_k1 - S_k0) > epsilon:
        k += 1
        S_k0 = S_k1
        S_k1 = ((-z)**k)/factorial(2*k+3)
        
    return S_k1

def delC(z,epsilon):
    
    k=2
    delC_0 = 1/factorial(2*k)
    delC_1 = delC_0 - k*((-z)**(k-1)/factorial(2*k+2)
    
    while abs(delC_1 - delC_0) > epsilon:
        k += 1
        delC_0 = delC_1
        delC_1 -= k*((-z)**(k-1))/factorial(2*k+2)
    
    return delC_1

def delS(z,epsilon):
    
    k=2
    delS_0 = 1/factorial(2*k+1)
    delS_1 = delS_0 - k*((-z)**(k-1)/factorial(2*k+3)
    
    while abs(delS_1 - delS_0) > epsilon:
        k += 1
        delS_0 = delS_1
        delS_1 -= k*((-z)**(k-1))/factorial(2*k+3)
    
    return delS_1

# Global Variables



# Main

