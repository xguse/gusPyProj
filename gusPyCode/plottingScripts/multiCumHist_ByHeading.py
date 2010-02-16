"""Produce an overlaid figure of transparent cumulative histograms
using the data headings as data labels.
"""

import optparse
from matplotlib import rc
from matplotlib import pylab as pl
rc('text', usetex=True)

print '\n\n'
usage = """python %prog inFile [options] """
parser = optparse.OptionParser(usage)

parser.add_option('-o', dest="outFile", default=None,
                  help="""Set outFile (default=Derive from inFile)""")
parser.add_option('-a', dest="alpha", default=0.6,type="float",
                  help="""Set alpha (default=%default)""")
parser.add_option('-t', dest="title", default="",type="string",
                  help="""Set title [TeX allowed] (default=%default)""")
parser.add_option('-n', dest="normed",action="store_true", default=False,
                  help="""Set to 'normed' hist (default=%default)""")
#parser.add_option('--bins',type='str', default="30",
                  #help="""Set Bins options (default=%default)""")
parser.add_option('--x_label',type='str', default='',
                  help="""Set x-axis label (default=%default)""")
parser.add_option('--y_label',type='str', default='',
                  help="""Set y-axis label (default=%default)""")
parser.add_option('--xMin',type='int', default=False,
                  help="""Set x-axis lower bound (default=%default)""")
parser.add_option('--xMax',type='int', default=False,
                  help="""Set x-axis upper bound (default=%default)""")
parser.add_option('--no-show', dest="no_show",action="store_true", default=False,
                  help="""Suppress figure displays (default=%default)""")


(opts, args) = parser.parse_args()

if len(args) == 0:
    parser.print_help()
    exit(0)
if len(args) != 1:
    print '**ERROR: Please supply only one inFile.\n'
    parser.print_help()
    exit(1)

if opts.outFile:
    outFile = opts.outFile
else:
    outFile = args[0]+'.pdf'
    

# --- Crunch File ---
data = map(lambda l: l.strip('\n').split('\t'), open(args[0], 'rU'))

heads = data.pop(0)


# --- Make Fig ---
fig = pl.figure()
fig.suptitle(r'%s' % (opts.title))
ax = fig.add_subplot(111)
ax.set_xlabel(r'%s' % (opts.x_label))
ax.set_ylabel(r'%s' % (opts.y_label))

for i in range(len(data[0])):
    tmpData = [x[i] for x in data]
    histData = []
    
    for datum in tmpData:
        if datum != '':
            histData.append(float(datum))
            
    #ax.hist(histData, bins=50, range=None, normed=opts.normed,
            #weights=None, cumulative=False, bottom=None,
            #histtype='stepfilled', align='mid', orientation='vertical',
            #rwidth=None, log=False, alpha=opts.alpha, label=heads[i])
    bins=[0.000,0.025,0.050,0.075,0.100,0.125,0.150,0.175,0.200,0.225,0.250,0.275,0.300,0.325,0.350,0.375,0.400,0.425,0.450,0.475,0.500,0.525,0.550,0.575,0.600]        
    ax.hist(histData, bins=bins, range=None, normed=opts.normed,
            weights=None, cumulative=True, bottom=None,
            histtype='stepfilled', align='mid', orientation='vertical',
            rwidth=None, log=False, alpha=opts.alpha, label=heads[i])
# --- Save/Show Fig ---
ax.legend(loc='upper left')
pl.savefig(outFile)
if not opts.no_show:
    pl.show()




