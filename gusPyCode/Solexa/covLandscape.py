"""Plot coverage landscapes.
"""
import matplotlib
matplotlib.use('Qt4Agg')
import sys
import optparse
from matplotlib import pylab as plb
print '\n\n\n'

    
usage = """python %prog inFile [options] """
parser = optparse.OptionParser(usage)
parser.add_option('-v',dest="plot_variations",action="store_true",default=False,
                  help="""Plots corresponding variation landscape (default=%default)""")
parser.add_option('-s', dest="plot_snpCalls", default=False,
                  help="""Plots SNP calls in provided filePath (default=%default)""")
parser.add_option('-o', dest="saveFig", action="store_true",default=False,
                  help="""Save figure (default=%default)""")
parser.add_option('-t',dest="title",type="string", default='',
                  help="""Quoted string: "the title text" """)
parser.add_option('--fdr', dest="fdrThresh", default=0.05,type="float",
                  help="""default=%default""")
parser.add_option('--no-show', dest="no_show",action="store_true", default=False,
                  help="""Suppress figure displays (default=%default)""")
parser.add_option('--xmin',type='int', default=False,
                  help="""Set x-axis lower bound (default=%default)""")
parser.add_option('--xmax',type='int', default=False,
                  help="""Set x-axis upper bound (default=%default)""")

(opts, args) = parser.parse_args()

if len(sys.argv) == 1:
    parser.print_help()
    exit(1)
if len(args) == 0:
    print "\n\nERROR: Please supply inFile."
    exit(1)
if opts.plot_snpCalls:
    print "\n\nNOTICE: I am using FDR=%s." % (opts.fdrThresh)
if opts.plot_snpCalls:
    opts.title+=' FDR=%s' % (opts.fdrThresh)

fig = plb.figure()
fig.suptitle(opts.title)
ax  = fig.add_subplot(1,1,1)

nucData = [[],[],[]]
inFile  = open(args[0],'rU')
print 'Reading file...'
for line in inFile:
    l = [int(x) for x in line.strip('\n').split('\t')]
    nucData[0].append(l[0])
    nucData[1].append(l[1])
    nucData[2].append(l[2])

# ------ Plot Coverage -------
print 'Plotting coverage...'
ax.plot(nucData[0],nucData[1],'b',mew=4,label='Coverage')
#ax.fill_between([1,2,3,4],[4,3,2,1],'b')


# ------ Plot Variations -------
if opts.plot_variations:
    print 'Plotting variation data...'
    ax.plot(nucData[0],nucData[2],'r',mew=4,label='Variation',alpha=0.5)
    #cAx.fill_between([1,2,3,4],[4,3,2,1],'b')


if opts.plot_snpCalls:
    print "Loading SNPs and FDRs..."
    snpFile = open(opts.plot_snpCalls,'rU')
    snpData = []
    for line in snpFile:
        l = [eval(x) for x in line.strip('\n').split('\t')]
        if l[5] <= opts.fdrThresh:
            snpData.append([l[0],l[2]])
        else: pass
# ------ Plot SNPs -------
    print "Plotting SNPs..."
    snpData.sort(key=lambda x: x[0])
    Xs = [x[0] for x in snpData]
    Ys = [y[1] for y in snpData]
    ax.bar(Xs,Ys,color='g',width=5,edgecolor='g',align='center',label='Putative SNPs Locations')


ax.legend(loc=0,shadow=True, fancybox=True)
plb.ylabel('Coverage')
plb.xlabel('Genomic Location')
if opts.xmin:
    ax.set_xbound(lower=opts.xmin)
if opts.xmax:
    ax.set_xbound(upper=opts.xmax)
if opts.saveFig:
    plb.savefig('%s.FIG.pdf'%(args[0]))

if not opts.no_show:
    #print opts.no_show
    plb.show()

    
    
