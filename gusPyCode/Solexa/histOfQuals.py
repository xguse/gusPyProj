"""Plot hists of quality scores.
"""

import matplotlib
matplotlib.use('Qt4Agg')
import sys
import optparse
from matplotlib import pylab as plb
print '\n\n\n'

    
usage = """python %prog [options]  inFile [xtraFile xtraFile]"""
parser = optparse.OptionParser(usage)

parser.add_option('-t',dest="title",type="string", default='', 
                  help="""Quoted string: "the title text" """)
parser.add_option('-o', dest="saveFig", action="store_true",default=False,
                  help="""Save figure (default=%default)""")
parser.add_option('--kwargs',dest="kwargs",type="string", default=None, 
                  help="""Keyword args appropriate for 'pylab.hist()' """)
parser.add_option('--no-show', dest="no_show",action="store_true", default=False,
                  help="""Suppress figure displays (default=%default)""")


(opts, args) = parser.parse_args()

if len(sys.argv) == 1:
    parser.print_help()
    exit(1)
if len(args) == 0:
    print "\n\nERROR: Please supply inFile."
    exit(1)


    
    
# ----- load data -----
matches = []
missAln = []

for f in args:
    for line in open(f):
        l = line.strip('\n').split('\t')
        if l[0] == 'm':
            matches.append(float(l[1]))
        elif l[0] == 'v':
            missAln.append(float(l[1]))

            
# ----- initiate figure -----
plb.rcParams['figure.figsize'] = 6, 11
fig = plb.figure()
fig.suptitle(opts.title)

# ----- plot aligned hist ----- 
ax  = fig.add_subplot(2,1,1)
ax.set_title('Aligned')

ax.hist(matches, bins=100, range=None, normed=1, cumulative=False,
        bottom=None, histtype='bar', align='mid',
        orientation='vertical', rwidth=None, log=1)


# ----- plot misAligned hist ----- 
fig.suptitle(opts.title)
ax2  = fig.add_subplot(2,1,2)
ax2.set_title('Miss-Aligned')

ax2.hist(missAln, bins=100, range=None, normed=1, cumulative=False,
        bottom=None, histtype='bar', align='mid',
        orientation='vertical', rwidth=None, log=1)

#ax.legend(loc=0,shadow=True, fancybox=True)
plb.ylabel('Base Calls')
plb.xlabel('Probability of Incorrect Call')

if opts.saveFig:
    plb.savefig('%s.HIST.pdf'%(args[0]))

if not opts.no_show:
    #print opts.no_show
    plb.show()

    
    

    
    

