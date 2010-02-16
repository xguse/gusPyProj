"""Calculate and plot histogram(s) from a table of data with headers:
#trascript_name	LS	LBBR	RS	RBA	CS	CB
AAEL000001-RA	279	404	481	331	1148	833
AAEL000002-RA	6	15	41	33	108	100
AAEL000003-RA	256	450	518	311	1523	998

Uses headers to determine which data to plot. Commented headers are required.

"""

import sys
import optparse
from matplotlib import rc
from matplotlib import pylab as pl
from gusPyCode.defs.JamesDefs import tbl2nmdTup,Bag
rc('text', usetex=True)

print '\n\n'

usage = """python %prog inFile [options] """
parser = optparse.OptionParser(usage)
parser.add_option('-i',dest="include",type='string',default=False,
                  help="""Quoted space delim'd list of headers to include (default=%default)""")
parser.add_option('-p',dest="prefix",type='string',default=False,
                  help="""Save figs using file-name prefix provided (default=inFile)""")
parser.add_option('--bins',dest="bins",type="string", default="30",
                  help="""Number of auto-bins to use, or quoted list of space delim'd left edges (default=%default)""")
#parser.add_option('--range',dest="range",type="string", default=None,
                  #help="""Quoted space delim'd "lowerBound upperBound" (default=%default)""")
parser.add_option('--xlabel',dest="xlabel",type="string", default='',
                  help="""Quoted string [TeX allowed]: "the xlabel text" (default=%default)""")
parser.add_option('--ylabel',dest="ylabel",type="string", default='',
                  help="""Quoted string [TeX allowed]: "the ylabel text" (default=%default)""")
parser.add_option('--cum',dest="cum",action='store_true', default=False,
                  help="""Plot as cumulative hist (default=%default)""")
parser.add_option('--log',dest="log",action='store_true', default=False,
                  help="""Plot as y-axis as log (default=%default)""")
parser.add_option('--no-show', dest="no_show",action="store_true", default=False,
                  help="""Suppress figure displays (default=%default)""")

(opts, args) = parser.parse_args()

if len(sys.argv) == 1:
    parser.print_help()
    exit(1)
if len(args) != 1:
    print "\n\nERROR: Please supply exactly one inFile."
    exit(1)
if not opts.prefix:
    opts.prefix = args[0]

opts.bins = opts.bins.split()
if len(opts.bins) == 1:
    opts.bins = int(opts.bins[0])
else:
    for i in range(len(opts.bins)):
        opts.bins[i] = int(opts.bins[i])

#if opts.range != None:
    #opts.range = opts.range.split()
    #if not (len(opts.range) == 2):
        #print "ERROR: range has len %s and not 2."
        #exit(1)
    #else:
        #opts.range[0],opts.range[1] = float(opts.range[0]),float(opts.range[1])
        #opts.range = (min(opts.range), max(opts.range))

  
# --- Parse inFile ---
table = tbl2nmdTup(args[0],'reads')
for i in range(len(table)):
    table[i] = Bag(table[i]._asdict())
if not opts.include:
    opts.include = table[0].keys()
else:
    opts.include = opts.include.split()

for i in range(len(table)):
    for column in opts.include: 
        try:
            table[i][column] = float(table[i][column])
        except ValueError:
            print 'ERROR: A column you asked me to include contains a data type other than "int" or "float"!'
            exit(1)
            
# ------- Build Figs -------
for column in opts.include:
    fig = pl.figure()
    ax  = fig.add_subplot(111)
    ax.set_title(r"{\Huge %s}" % (column))
    ax.set_xlabel(r"{\Large %s}" % (opts.xlabel))
    ax.set_ylabel(r"{\Large %s}" % (opts.ylabel))
    
    data  = [x[column] for x in table]
    # ---- Extract Zeros -----
    zeros  = []
    nZeros = []
    for i in range(len(data)):
        if data[i] == 0:
            zeros.append(data[i])
        else:
            nZeros.append(data[i])
    data = nZeros
    
    n, bins, patches  =ax.hist(data, bins=opts.bins, range=[1,max(data)], normed=False,
                               cumulative=opts.cum, bottom=None, histtype='bar',
                               align='mid', orientation='vertical', rwidth=None,
                               log=opts.log,)
    ax.set_xscale('log')
    ax.text(0.8, 0.95, 
            'Zero = %s' % (len(zeros)),
            bbox=dict(facecolor='grey', alpha=0.5),
            transform = ax.transAxes)
    
    
    pl.savefig('%s.%s.hist.pdf' % (opts.prefix,column))
    if not opts.no_show:
        pl.show()
    del(fig,ax)
    