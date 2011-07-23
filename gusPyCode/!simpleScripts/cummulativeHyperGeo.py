from math import log10
import sys
import optparse
import csv
import collections


def benjHochFDR(pVals,pValColumn=1,FDR=0.05):
    """
    pVals      = 2D list(hypothesis,p-value) hypothesis could = geneName tested for enrichment
    pValColumn = integer of column index containing the p-value.
    !*! FDR    = threshold above which FDR is unacceptable [not yet implemented]!*! 
    
    Returns, for all *acceptable q-values: hypothesis,origPval,adjustedPval 
    *NOTE:  For now, returns _ALL_ items passed to it with no filtering at the moment.
    """
    assert type(pValColumn) == type(1),\
           "ERROR: pValColumn must be int type!"
    # Sort pVals from highest to lowest after converting them all to floats.
    for i in range(len(pVals)):
        pVals[i][pValColumn] = float(pVals[i][pValColumn])
    pVals.sort(key=lambda x: x[pValColumn])
    pVals.reverse()
    
    n = len(pVals)
    
    lastPval = pVals[0][pValColumn]
    for i in range(len(pVals)):
        p    = pVals[i][pValColumn]
        adj  = (float(n)/(n-i))
        adjP = p*adj
        miN  = min(adjP,lastPval)
        pVals[i].append(miN)
        lastPval = pVals[i][-1]
    
    pVals.reverse()
    return pVals

def tableFile2namedTuple(tablePath,sep='\t'):
    """Returns namedTuple from table file using first row fields as col headers."""
    #import collections
    #import csv
    
    reader  = csv.reader(open(tablePath), delimiter=sep)
    headers = reader.next()
    Table   = collections.namedtuple('Table', ', '.join(headers))
    data    = map(Table._make, reader)
    return data

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
    return 10**(log10(choose(n,i))+log10(choose(m,N-i))-log10(choose(n+m,N)))





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

    upperLim = min(N,n)
    
    for x in range(i,upperLim+1):
        cumPVal = cumPVal + hypergeoP(n,x,m,N)

    return cumPVal


def updateOutDict(outDict,rowData):
    """Create and/or add rowData to dict key row.CorrectionGroup
    which is at rowData[0].
    """
    try:
        outDict[rowData[0]].append(rowData)
    except KeyError:
        outDict[rowData[0]] = []
        outDict[rowData[0]].append(rowData)
        
def getQvals(outDict,pValColumn):
    """Call benjHochFDR() on all value lists in outDict,
    doing the work in place.
    """
    for group in outDict:
        benjHochFDR(outDict[group],pValColumn)

if __name__ == "__main__":

    
    #+++++++++++ File Parseing Etc +++++++++++
    columnNames = ['CorrectionGroup',
                   'SampleName',
                   'PositivesInPopulation',
                   'PositivesInSample',
                   'NegativesInPopulation',
                   'SampleSize']
    extraHelp = """\n\nNote: INFILE must contain the following EXACT column titles
    in the first row (column order does not matter):
    %s\n\n""" % ('", "'.join(columnNames))
    
    usage = """python %prog -i inFile -o outFile"""
    parser = optparse.OptionParser(usage=usage)
    
    parser.add_option('-i', dest="infile", type='string',default=None,
                      help="""Path to inFile. (default=%default)""")
    parser.add_option('-o', dest="outfile", type='string',default=None,
                      help="""File to write results to. (default=%default)""")
    parser.add_option('--delim', dest="delim", type='str',default='\t',
                      help="""What character is used to delimit the inFile? (default="\\t")""")
    
    if ('-h' in sys.argv) or ('--help' in sys.argv):
        print extraHelp
        
    (opts, args) = parser.parse_args()
    
    
    if len(sys.argv) == 1:
        print extraHelp
        parser.print_help()
        exit(0)
    if (opts.infile == None) or (opts.outfile == None):
        raise Exception("ERROR: you must provide both INFILE and OUTFILE.")
    
    # Load table file and check for correct header names
    table = tableFile2namedTuple(tablePath=opts.infile,sep=opts.delim)
    tableFields = table[0]._fields
    for f in tableFields:
        if f not in columnNames:
            raise Exception("ERROR: column name '%s' does not match %s" % (f, columnNames))
    if len(tableFields) != len(columnNames):
        raise Exception("ERROR: your file has %s columns while it should have exactly %s columns." %
                        (len(tableFields),len(columnNames)))
    
    # Initialize output list
    outputHeaders = columnNames + ['UnCorrectedPvalues','BH:qVals']
    
    output = {}
    # run the analysis:
    for row in table:
        n=int(row.PositivesInPopulation)
        i=int(row.PositivesInSample)
        m=int(row.NegativesInPopulation)
        N=int(row.SampleSize)
        if ((n > 0) and (i > 0)):
            pVal = cumHypergeoP(n,i,m,N)
            
            updateOutDict(output,[row.CorrectionGroup,
                                  row.SampleName,
                                  row.PositivesInPopulation,
                                  row.PositivesInSample,
                                  row.NegativesInPopulation,
                                  row.SampleSize,
                                  str(pVal)])
        
    
    getQvals(output,-1)
    
    # write the output
    outFile = open(opts.outfile,'w')
    outFile.write('%s\n' % ('\t'.join(outputHeaders)))
    output = output.values()
    outLines = []
    for group in output:
        outLines.extend(group)
    outLines.sort(key=lambda x: x[0])
    
    for line in outLines:
        outFile.write('%s\n' % ('\t'.join([str(x) for x in line])))
    outFile.close()