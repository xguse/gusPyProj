import optparse
import sys



def updateCounts(term,countDict):
    """"""
    try:
        countDict[term] += 1
    except KeyError:
        countDict[term] = 1

if __name__ == "__main__":

    
    #+++++++++++ File Parseing Etc +++++++++++
    epilog = """"""
    
    usage = """python %prog [options] -i infFile -c columnNumber"""
    parser = optparse.OptionParser(usage=usage, epilog=epilog)
    
    parser.add_option('-i', dest="infile", type='string',default=None,
                      help="""Path to inFile. (default=%default)""")
    parser.add_option('-c', dest="column", type='int',default=1,
                      help="""Which column to count. (default=%default)""")
    parser.add_option('--delim', dest="delim", type='str',default='\t',
                      help="""What character is used to delimit the inFile? (default="\\t")""")
    
    
    (opts, args) = parser.parse_args()
    
    
    if len(sys.argv) == 1:
        parser.print_help()
        exit(0)
    
    inFile = open(opts.infile,'rU')
    
    # do the counting
    countDict = {}
    for line in inFile:
        line = line.strip('\n').split(opts.delim)
        updateCounts(line[opts.column - 1],countDict)
        
    # sort terms
    terms = countDict.keys()
    terms.sort()
    
    # print results
    for t in terms:
        print "%s\t%s" % (t,countDict[t])
        