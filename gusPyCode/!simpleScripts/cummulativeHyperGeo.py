import sys
import optparse
import csv
import collections

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
    return (choose(n,i)*choose(m,N-i))/float(choose(n+m,N))




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


if __name__ == "__main__":

    
    #+++++++++++ File Parseing Etc +++++++++++
    columnNames = ['SampleName',
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
    output = [columnNames + ['UnCorrectedPvalues']]
    
    # run the analysis:
    for row in table:
        pVal = cumHypergeoP(n=int(row.PositivesInPopulation),
                            i=int(row.PositivesInSample),
                            m=int(row.NegativesInPopulation),
                            N=int(row.SampleSize))
        output.append([row.SampleName,
                       row.PositivesInPopulation,
                       row.PositivesInSample,
                       row.NegativesInPopulation,
                       row.SampleSize,
                       str(pVal)])
    
    # write the output
    outFile = open(opts.outfile,'w')
    for line in output:
        outFile.write('%s\n' % ('\t'.join(line)))
    outFile.close()