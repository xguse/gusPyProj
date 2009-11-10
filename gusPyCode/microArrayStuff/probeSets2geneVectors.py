print "loading modules..."

from scipy.stats.stats import pearsonr
#from scipy.spatial.distance import pdist
from defs_microArray import geneVectors
from gusPyCode.defs.JamesDefs import removeCommentLines


def pairwiseMetric



def main(probSetVectors, probSet2geneList, metric='pearson', cutOff='default'):
    """Takes a list of prob sets and associated expression vectors, 
    a list of probe set-to-gene associations, and user options.  Returns
    a list of genes with associated expression vectors calculated from
    the set of probeSets associated with the gene using the median. 
    User options: 
    distance metric (pearson -or- euclidean)
    """
    print 'starting "main"...'
    
    # Initiate instance of geneVectors class for probsets
    probSets = geneVectors(removeCommentLines(probSetVectors, '#'))
    
    # Format probSet2geneList as list of lists
    probSet2geneList = map(lambda l: l.split('\t'), probSet2geneList)
        
    # Initiate dictionary of gene names to accept probset ids
    genes    = {}
    for i in probSet2geneList:
        genes[i[1]] = []
        
    # Group probsets by gene recording and eliminating those that
    # match more than one gene.
    wantonPSs = []
    for psID in probSets.vectors.keys():
        g = []
        for association in probSet2geneList:
            if association[0] == psID:
                g.append(association[1])
        
        
        if len(g) > 1:  # collect non-specific probes and the genes they are associated with for report
            wantonPSs.append([psID,g[:]])
        elif len(g) == 1:
            genes[g[0]].append(psID)
        elif len(g) == 0: 
            print 'warn: %s matched no genes.' % (psID) 
            
    # If more than one probSet, measure and attached metric
    for gene in genes.keys():
        

            
    pass







if __name__ == '__main__':
    import sys
    
    # Parse args
    args = sys.argv[1:]
    assert len(args) >= 4, 'usage: %s probSetVectorFile probSet2geneListFile metric<"pearson"> cutOff<"default"> [optionalOutFile]' % (sys.argv[0])
    
    # Load files
    probSetVectors   = map(lambda line: line.strip(), open(args[0], 'rU').readlines())
    probSet2geneList = map(lambda line: line.strip(), open(args[1], 'rU').readlines())
    
    # Feed lists to main
    results = main(probSetVectors, probSet2geneList, args[2], args[3])
    
    # Route the output based on whether an outfile was supplied
    if args[3]: ## write to file
        oFile = open(args[3], 'w')
        for i in results:
            oFile.write('%s\n' % (i))
        oFile.close()
    else: ## print to screen
        for i in results:
            print i
    
    None

