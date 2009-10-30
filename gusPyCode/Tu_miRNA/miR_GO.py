from gusPyCode.defs.statsDefs import hypergeoP



def goEnrichment(miR_targets,GOgenes,popSize):
    """
    goEnrichment(miR_targets,GOgenes,popSize):
    miR_targets = set(predictedGenes)
    GOgenes     = set(genesInGOgroup)
    popSize     = int(numberOfGenesConsideredAsPopulation)
    
    Returns cumulative hypergeometric enrichment P-value.
    """
    #n = # of positives in population
    #i = # of positives in sample
    #m = # of negatives in population
    #N = sample size
    #P(x=i) = (choose(n,i)choose(m,N-i))/choose(n+m,N)
    #For more details -> http://mathworld.wolfram.com/HypergeometricDistribution.html
    n = len(GOgenes)
    i = len(GOgenes.intersection(miR_targets))
    m = popSize-len(GOgenes)
    N = len(miR_targets)
    
    return sum([hypergeoP(n,x,m,N) for x in range(i,N+1)])
    