from time import time
from time import ctime
import optparse
import sys
import os
import scipy.stats as stats
import matplotlib as mpl
if os.environ['USER'] == 'biggus':
    mpl.use('TkAgg')
try:
    import gusPyCode
except:
    sys.path.extend(['/Users/biggus/lib/python/site-packages/',
                     '/Library/Python/2.5/site-packages/',
                     '/Library/Frameworks/Python.framework/Versions/2.5/lib/python2.5/site-packages/'])
from matplotlib import pylab as pl
from gusPyCode.defs.JamesDefs import Bag
from gusPyCode.defs.bioDefs import ParseFastQ, ParseBowtieBed


 

#+++++++++ Definitions ++++++++++ 

def parseVecfiles(pathsList):
    """Takes pathsList and returns a list of parsed vectorFile info objects.
    vectorInfoObj format = tuple(t(header0,header1),t(vector0),t(vector1))"""
    
    rList  = []
    
    for path in pathsList:
        header = []
        vec0   = []
        vec1   = []
        
        vecFile = open(path,'rU')
        
        header.extend(vecFile.readline().strip('\n').split('\t'))
        #dBugCnt = 0
        while 1: #dBugCnt < 10000:
            #dBugCnt+=1
            line = vecFile.readline()
            if not line:
                break
            line = line.strip('\n').split('\t')
            vec0.append(int(line[0]))
            vec1.append(int(line[1]))
        rList.append(tuple([tuple(header),tuple(vec0),tuple(vec1)]))
    
    return tuple(rList)
        
            
            
        
    


if __name__ == '__main__':
    print '\n\n\n'
    
    #+++++++++ Parse Command Line ++++++++++
        
    usage = """python %prog [options] vectorFile0 [vectorFile1 vectorFile2 ...]"""
    parser = optparse.OptionParser(usage)

    parser.add_option('--show',dest="show",action="store_true",default=False,
                      help="""Show plot(s) in window. (default=%default)""")


    
    
    (opts, args) = parser.parse_args()
    
    if len(sys.argv) == 1:
        parser.print_help()
        exit(1)
    if len(args) < 1:
        parser.print_help()
        print "\n\n** ERROR: Please supply at least one vectorFile. **"
        exit(1)

        

    # ++++++++++ Open Files And Get Variables Set up ++++++++
    print 'Parsing vector files...' 
    vectorFiles   = parseVecfiles(args)


    # ++++++++++ plot each vector file ++++++++
    for vecFile in vectorFiles:
        
        # ++++++++++ Calculate The Correlations ++++++++ 
        
        print 'Calculating Pearson... %s' % (ctime())
        pearson  = stats.pearsonr(vecFile[1],vecFile[2])
        print pearson
    
        
        # ++++++++++ Plot the scatter plot ++++++++ 
        print "Drawing..."
        fig = pl.figure()
        ax  = fig.add_subplot(111)
        ax.set_xscale('log')
        ax.set_yscale('log')
        
        ax.scatter(vecFile[1],vecFile[2], s=15, c='b', marker='o', alpha=0.5)
        ax.set_xlabel(vecFile[0][0])
        ax.set_ylabel(vecFile[0][1])

        
        m,b  = pl.polyfit(vecFile[1],vecFile[2],1)
        bfYs = pl.polyval([m,b], [1,max(vecFile[1])])
        
        ax.plot([1,max(vecFile[1])],bfYs,'r-')
        
        pl.text(0.01,0.99,'Pearson: %.4f, %s\nBest Fit: y=%.3f*x+%.3f' % (pearson[0],pearson[1],m,b),
                bbox=dict(facecolor='#87AACD', alpha=1),
                horizontalalignment='left',
                verticalalignment='top',
                transform = ax.transAxes)
        
        

        pl.savefig('%s_vs_%s.png' % (vecFile[0][0],vecFile[0][1]))
        print 'Show?  %s' % (opts.show)
        if opts.show:
            pl.show()
    