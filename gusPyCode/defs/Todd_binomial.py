#!/usr/bin/python

# Binomial coefficients. Found on interwebz.
def binc(n, k):
    """binc(n, k): Computes n choose k."""
    if (k > n): return 0
    if (k < 0): return 0
    if (k > int(n/2)):
        k = n - k
    rv = 1
    for j in range(0, k):
        rv *= n - j
        rv /= j + 1
    return rv

def pk(k,n,ko,no):
    p = (ko/float(no));
    return binc(n,k) * (p**k) * ((1-p)**(n-k));

# ki is number of appearances in G
# n is size of G
# ko is number of appearances in background
# no is size of background
def pvalue(ki,n,ko,no) :
    tot = 0;
    for i in range(ki,n+1):
        tot += pk(i,n,ko,no);
    return tot;

ki = 5
n = 10
ko = 200
no = 3000

print binc(n, ki)
print pk(ki, n, ko, no)
print pvalue(ki, n, ko, no)