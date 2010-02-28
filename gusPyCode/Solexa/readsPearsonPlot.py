from time import time
from time import ctime
import optparse
import sys
import scipy.stats as stats
import matplotlib as mpl
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

supportedFileTypes = {"bowtie_bed":ParseBowtieBed,
                      "fastq": ParseFastQ}
def setParser(filePath,fileType):
    """Return the correct parser initialized and ready for use."""
    assert fileType in supportedFileTypes,\
           "** ERROR: %s not supported file type: %s **" % (fileType,supportedFileTypes.keys())
    return supportedFileTypes[fileType](filePath)

    
def addRead(readSeq,readDict,whichFile):
    """Increment count for read in correct file count. If not in dict,
    initialize entry and increment correct file's count."""
    
    if readSeq in readDict:
        readDict[readSeq][whichFile] += 1
    else:
        readDict[readSeq] = [0,0]
        readDict[readSeq][whichFile] += 1


if __name__ == '__main__':
    print '\n\n\n'
    
    #+++++++++ Parse Command Line ++++++++++
        
    usage = """python %prog -o outName readsFile1 readsFile2"""
    parser = optparse.OptionParser(usage)
    parser.add_option('-o',dest="out_name",type="string",default=False,
                      help="""<required> Name to give result files.""")
    parser.add_option('-f',dest="file_type",type="string",default='fastq',
                      help="""File format of readsFiles.""")
    parser.add_option('--show',dest="show",action="store_true",default=False,
                      help="""Show plot in window.""")


    
    
    (opts, args) = parser.parse_args()
    
    if len(sys.argv) == 1:
        parser.print_help()
        exit(1)
    if len(args) != 2:
        parser.print_help()
        print "\n\n** ERROR: Please supply exactly two readsFiles. **"
        exit(1)
    if not opts.out_name:
        parser.print_help()
        print "\n\n** ERROR: You must supply an out name. **"
        

    # ++++++++++ Open Files And Get Variables Set up ++++++++
    readsFile0   = args[0]
    readsFile1   = args[1]
    # -- set up correct parsers --
    parser0 = setParser(readsFile0,opts.file_type)
    parser1 = setParser(readsFile1,opts.file_type)
    
    readDict     = {} ## to collect read counts
    
    # ++++++++++ Execute the counting ++++++++
    # --- readsFile0 ---
    print 'Counting the first file...'
    t1_0 = time()
    whichFile = 0
    while 1:
        readSeq = parser0.getNextReadSeq()
        if readSeq:
            addRead(readSeq,readDict,whichFile)
        else: break
    t2_0 = time()
    print 'Counting took %s min.' % ((t2_0-t1_0)/60)
    # --- readsFile1 ---
    print 'Counting the second file...'
    t1_1 = time()
    whichFile = 1
    while 1:
        readSeq = parser1.getNextReadSeq()
        if readSeq:
            addRead(readSeq,readDict,whichFile)
        else: break
    t2_1 = time()
    print 'Counting took %s min.' % ((t2_1-t1_1)/60)
        
    # ++++++++++ Calculate The Correlations ++++++++ 
    print 'Getting vectors... %s' % (ctime())
    vector0 = tuple([x[0] for x in readDict.values()])
    vector1 = tuple([x[0] for x in readDict.values()])
    print ctime()
    
    print 'Deleting readDict...'
    t1_del = time()
    del(readDict)
    t2_del = time()
    print '%s min.' % ((t2_del-t1_del)/60)
    
    print 'Calculating Pearson... %s' % (ctime())
    pearson  = stats.pearsonr(vector0,vector1)
    print 'Calculating Spearman... %s' % (ctime())
    spearman = stats.spearmanr(vector0,vector1)
    
    # ++++++++++ Plot the scatter plot ++++++++ 
    print "Drawing..."
    fig = pl.figure()
    ax  = fig.add_subplot(111)
    
    ax.plot(vector0,vector1, 'bo')
    ax.set_xlabel(parser0._file.name)
    ax.set_ylabel(parser1._file.name)
    pl.text(0,1,'Pearson(corr,p-val[2t]): (%s,%s)\nSpearman(corr,p-val[2t]): (%s,%s)' % (pearson[0],pearson[1],spearman[0],spearman[1]),
            bbox=dict(facecolor='grey', alpha=1),
            horizontalalignment='left',
            verticalalignment='top',
            transform = ax.transAxes)
    pl.savefig(opts.out_name)
    print 'Show?  %s' % (opts.no_show)
    if opts.show:
        pl.show()
    