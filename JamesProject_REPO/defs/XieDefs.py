#!/usr/bin/python
"""
mkStats.py

Matt Kayala
CS 284A 

Stats Utilities for Assignment 1


A function to calculate p-value in a hypergeometric distr.

Have to make some approximation by working with logs.

"""

import os, sys;
from math import *;


def choose(n, k):
    """
    Function for (n choose k)
    """
    return factorial(n)/factorial(k)/factorial(n-k);
    
    
def factorial(n):
    """
    Function for factorial.
    """
    if n == 0: return 1;
    fact = 1;
    for i in xrange(2, n+1):
        fact *= i;
    return fact;
    
def hyperGeom(N, n, K, k):
    """
    Function for calcing hypergeometric dist prob mass function
    
    Uses log so sacrifices a little bit of accuracy.
    """
    hyp=0;
    try:
        hyp = log(choose(K, k)) + log(choose(N-K, n-k)) - log(choose(N, n));
    except ValueError, inst:
        print 'K', K;
        print 'N', N;
        print 'k', k;
        print 'n',  n;
        raise inst;
        
    return exp(hyp);
    
    
def pValHypGeom(N, n, K, k):
    """
    Simple function for calculating \sum p(i) for i from k to n
    """
    pVal = 0;
    limit = min(n, K);
    
    for i in xrange(k, limit +1):
        #print(i);
        pVal += hyperGeom(N, n, K, i);
        
    return pVal;


k = 5
n = 10
K = 25
N = 120

pValVal = pValHypGeom(N,n,K,k)

print pValVal

