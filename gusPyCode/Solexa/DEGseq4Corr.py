from time import time
from time import ctime
import optparse
import sys
import os
import scipy.stats as stats
import matplotlib as mpl
mpl.use('TkAgg')
from matplotlib import pylab as pl
from gusPyCode.defs.JamesDefs import Bag
from gusPyCode.defs.bioDefs import ParseDEGseqOut


 

#+++++++++ Definitions ++++++++++ 
def toInt(string):
    if string == 'NA':
        return 0
    else:
        return int(string)


if __name__ == '__main__':
    print '\n\n\n'
    
    #+++++++++ Parse Command Line ++++++++++
        
    usage = """python %prog [options] DEGseqFile1 [DEGseqFile2 DEGseqFile3 ...]"""
    parser = optparse.OptionParser(usage)
    parser.add_option('--show',dest="show",action="store_true",default=False,
                      help="""Show plot(s) in window. (default=%default)""")
    


    
    
    (opts, args) = parser.parse_args()
    
    if len(sys.argv) == 1:
        parser.print_help()
        exit(1)
    if (len(args) < 1):
        parser.print_help()
        print "\n\n** ERROR: Please supply at least one DEGseqFile **"
        exit(1)


    # ++++++++++ Open Files And Get Variables Set up ++++++++
    for f in args:

        # -- set up file parser --
        DEGseqParser = ParseDEGseqOut(f)
    
        # -- get toal reads in each data-set --
        tots1 = 0
        tots2 = 0
        for line in DEGseqParser.data:
            tots1 += toInt(line[1])
            tots2 += toInt(line[2])
        
        
        # ++++++++++ Calculate The Correlations ++++++++

    
        
        print 'Getting vectors...' 
        vector0 = tuple([float(toInt(x[1]))/tots1 for x in DEGseqParser.data])
        vector1 = tuple([float(toInt(x[2]))/tots2 for x in DEGseqParser.data])

        

        
        print 'Calculating Pearson...'
        pearson  = stats.pearsonr(vector0,vector1)
        print 'max(xVect,yVect): %s,%s' % (max(vector0),max(vector1))
        print pearson
        
        # ++++++++++ Plot the scatter plot ++++++++ 
        print "Drawing..."
        fig = pl.figure()
        ax  = fig.add_subplot(111)
        #ax.set_xscale('log')
        #ax.set_yscale('log')
        
        ax.scatter(vector0,vector1, s=10, c='b', marker='o', alpha=0.6)
        
        # -- set axis labels --
        xLab = DEGseqParser._file.name.split('_')[0]
        yLab = DEGseqParser._file.name.split('_')[2]
        ax.set_xlabel(xLab)
        ax.set_ylabel(yLab)

        
        m,b  = pl.polyfit(vector0,vector1,1)
        min_xVec = min(vector0)
        max_xVec = max(vector0)
        bfYs = pl.polyval([m,b], [min_xVec,max_xVec])
        
        ax.plot([min_xVec,max_xVec],bfYs,'r-')
        
        pl.text(0.01,0.99,'Pearson: %.4f, %s\nBest Fit: y=%.3f*x+%.3f' % (pearson[0],pearson[1],m,b),
                bbox=dict(facecolor='#87AACD', alpha=1),
                horizontalalignment='left',
                verticalalignment='top',
                transform = ax.transAxes)
        
        

        pl.savefig('%s_vs_%s.indvDEG.png' % (xLab,yLab))
        print 'Show?  %s' % (opts.show)
        if opts.show:
            pl.show()       
    
        
        