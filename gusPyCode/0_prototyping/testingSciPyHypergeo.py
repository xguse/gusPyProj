from scipy.stats import hypergeom

"""
Hypergeometric distribution

    Models drawing objects from a bin. M is total number of objects, n is total number of Type I objects.
    RV counts number of Type I objects in N drawn without replacement from population.

    hypergeom.pmf(k, M, n, N) = choose(n,k)*choose(M-n,N-k)/choose(M,N) for N - (M-n) <= k <= min(m,N)
"""

k = 10         # [wiki=k] Targets found in drawn sample
M = 1000       # [wiki=N] Total objects in population
n = 30         # [wiki=m] Total targets in population
N = 50         # [wiki=n] Size of sample drawn

#dhyper(sampHits_x=10,totHits_m=30,totMiss_n=1000-30,)

pmf = hypergeom.pmf(k, M, n, N)
prb = hypergeom.cdf(k, M, n, N)

None