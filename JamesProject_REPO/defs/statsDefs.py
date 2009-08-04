from decimal import Decimal
# see bottom for conditional import of "bestChoose"



def setChoose():
    try:
	from gmpy import comb
    
    except ImportError:
	return binComb
    else:
	return comb


def binComb(n, k):
    """
    binComb(n, k): Computes n choose k. Defines the number of k objects that can be chosen from 
    n objects.
    """
    if (k > n): return 0
    if (k < 0): return 0
    if (k > int(n/2)):
	    k = n - k
    rv = 1
    for j in range(0, k):
	    rv *= n - j
	    rv /= j + 1
    return rv


def binom(n, m):
    """
    Calulates n-choose-k.
    
    def binom(n, m) is from "Data Structures and Algorithms
    with Object-Oriented Design Patterns in Python"
    by Bruno R. Preiss.

    Copyright (c) 2003 by Bruno R. Preiss, P.Eng. All rights reserved.

    http://www.brpreiss.com/books/opus7/programs/pgm14_10.txt
    """
    b = [0] * (n + 1)
    b[0] = 1
    for i in xrange(1, n + 1):
        b[i] = 1
        j = i - 1
        while j > 0:
            b[j] += b[j - 1]
            j -= 1
    return b[m]

def choose(n, k):
    if (k > n): return 0
    if (k < 0): return 0
    ntok = 1
    for t in xrange(min(k, n-k)):
        ntok = ntok*(n-t)//(t+1)
    return ntok

def hypergeoP(n,i,m,N):
    """
    Calculates the non-cumulative hypergeometric p-value for variables:
    n = # of positives in population
    i = # of positives in sample
    m = # of negatives in population
    N = sample size
    
    P(x=i) = (choose(n,i)choose(m,N-i))/choose(n+m,N)
    
    For more details -> http://mathworld.wolfram.com/HypergeometricDistribution.html
    """
    return (bestChoose(n,i)*bestChoose(m,N-i))/float(bestChoose(n+m,N))


    

def cumHypergeoP(n,i,m,N):
    """
    Calculates the cumulative hypergeometric p-value for variables:
    n = # of positives in population
    i = # of positives in sample
    m = # of negatives in population
    N = sample size
    
    P(i) = sum([as i->N] (choose(n,i)choose(m,N-i))/choose(n+m,N))
    
    For more details -> http://mathworld.wolfram.com/HypergeometricDistribution.html
    """
    
    cumPVal = 0
    
    for x in range(i,N+1):
	cumPVal = cumPVal + hypergeoP(n,x,m,N)
    
    return cumPVal





# determine whether we can use gmpy and set "bestChoose" to that or the standby.
try:
    from gmpy import comb as bestChoose
except ImportError:
    bestChoose = binComb
