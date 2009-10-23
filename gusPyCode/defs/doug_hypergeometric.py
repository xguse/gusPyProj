#!/usr/bin/python


# Binomial coefficients. Found on interwebz.
"""binc(n, k): Computes n choose k. Defines the number of k objects that can be chosen from 
n objects. For example, k could be "ki" and n "no" as specified below"""
def binc(n, k):
	if (k > n): return 0
	if (k < 0): return 0
	if (k > int(n/2)):
		k = n - k
	rv = 1
	for j in range(0, k):
		rv *= n - j
		rv /= j + 1
	return rv

def hypergeo(ki,n,ko,no):
	from decimal import Decimal
	""" number of ways of choosing k sequences with a motif from the total
   number ko sequences that also have the motif"""
	with_motif = binc(ko, ki)
	""" number of ways of choosing n-k sequences without a motif from
   the total number no-ko sequences also without the motif """
	without_motif = binc((no-ko),(n-ki))
	""" the number of ways of choosing n sequences from the total number
   N sequences """
	total_sequences = binc(no, n)
	returnThis = (with_motif*without_motif)/Decimal(total_sequences)
	#print returnThis
	return returnThis


# ki is number of appearances in "promoters" of a specific cluster
# n is number of "promoters" in the specific cluster
# ko is number of appearances in all "promoters"
# no is total number of "promoters"

#ki=5
#n=10
#ko= 25
#no= 50


def hyperGeoPvalue(no,n,ko,ki):
	tot = 0
	for i in range(ki,n+1):
		hypGeo = hypergeo(i,n,ko,no)
		tot = tot + hypGeo
		#print tot
	return tot








#print "the discrete hypergeometric probability:", hypergeo(ki, n, ko, no)
#pvalVal = pvalue(ki, n, ko, no)
#print "the cumulative probability:", pvalVal
#print range(ki, n+1)