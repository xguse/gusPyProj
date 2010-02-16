"""
Read in a table of reported FDR_Cas and produce a boxPlot representing FDR_Ca for each seed-type.
"""
import sys
import optparse
from matplotlib import rc
from matplotlib import pylab as pl
from collections import namedtuple
from gusPyCode.defs.JamesDefs import Bag
from gusPyCode.defs.mpl_custom import violin_plot
rc('text', usetex=True)

print '\n\n'

usage = """python %prog inFile [options] """
parser = optparse.OptionParser(usage)
parser.add_option('-b',dest="box_plot",action="store_true",default=False,
                  help="""Superimpose box-plot (default=%default)""")
parser.add_option('-o',dest="out_file",type='string',default=False,
                  help="""Save fig using file name provided (default=%default)""")
parser.add_option('-t',dest="title",type="string", default='',
                  help="""Quoted string: "the title text" """)
parser.add_option('--no-show', dest="no_show",action="store_true", default=False,
                  help="""Suppress figure displays (default=%default)""")

(opts, args) = parser.parse_args()

if len(sys.argv) == 1:
    parser.print_help()
    exit(1)
if len(args) != 1:
    print "\n\nERROR: Please supply exactly one inFile."
    exit(1)

# ---- Parse InFile ----
data = map(lambda l: l.strip('\n').split('\t'), open(args[0]).readlines())

if data[0][0].startswith('#'):
    data[0][0] = data[0][0][1:]
    header = data.pop(0)
else:
    print 'ERROR:  The first line of inFile must be the header and begin with a "#".'
    exit(1)
seedTypeData = namedtuple('seedTypeData',' '.join(header))
tabledata = [seedTypeData._make(x) for x in data]

# ---- Gather Data ----
bag = Bag()
for item in header:
    bag[item]=[]

for entry in tabledata:
    for sType in entry._fields:
        try:
            bag[sType].append(float(entry._asdict()[sType]))
        except ValueError:
            pass


# --- Build Fig ---

fig = pl.figure()
ax = fig.add_subplot(111)
if opts.title:
    ax.set_title(r'{\Huge %s}' % (opts.title))
ax.set_xlabel(r'{\huge Seed-Type}')
ax.set_ylabel(r'{\huge FDR$_{C_{a}}$}')
ax.set_xticklabels([x.replace('_','\_') for x in header])

vectors = []
for x in header:
    if bag[x]:
        vectors.append(bag[x])
    else:
        vectors.append(0)
    
violin_plot(ax,vectors,range(len(header)), bp=opts.box_plot)

if opts.out_file:
    pl.savefig(opts.out_file)
if not opts.no_show:
    pl.show()




print "Done"

        
    
