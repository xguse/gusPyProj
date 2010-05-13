import sys
import optparse
from gusPyCode.defs.JamesDefs import tableFile2namedTuple

strandReps = {'+':'+',
              '-':'-',
              '1':'+',
              '-1':'-',}

def groupFeatBlks(featBlkTable,opts):
    """Returns Dict with key,val mappings == feature_ID,[rowsFrom_featBlkTable]"""
    featureDict = {}
    for row in featBlkTable:
        if row.__getattribute__(opts.featName) in featureDict:
            featureDict[row.__getattribute__(opts.featName)].append(row)
        else:
            featureDict[row.__getattribute__(opts.featName)] = []
            featureDict[row.__getattribute__(opts.featName)].append(row)
    return featureDict



def printBEDline(listOfrowsByFeature,opts):
    """Takes a list of block info rows grouped by a single feature.
    Prints the BED line calculated from this information."""
    
    # +++++ func specific Defs +++++
    def getBlockSizes(feat):
        blkSzList = []
        for blk in feat:
            blkSzList.append(str(int(blk.__getattribute__(opts.blkChmEnd))-int(blk.__getattribute__(opts.blkChmStrt))+1))
        return ','.join(blkSzList)
    
    def getBlockStarts(feat,chrmStart):
        blkStrtList = []
        for blk in feat:
            blkStrtList.append(str(int(blk.__getattribute__(opts.blkChmStrt))-1-int(chrmStart)))
        return ','.join(blkStrtList)
    
    # +++ Do some data validation +++
    # --- Check Chrom info ---
    chrms = set([x.__getattribute__(opts.chrm) for x in listOfrowsByFeature])
    if len(chrms) != 1:
        print >> sys.stderr, 'ERROR: the %s feature has blocks on more than one chromosome/contig %s (%s column).  SKIPPING' % \
              (listOfrowsByFeature[0].__getattribute__(opts.featName),
               list(chrms),
               opts.chrm)
        return None
        #raise Exception, \
              #'ERROR: the %s feature seems to have blocks on more than one chromosome/contig\n%s (%s column)' % \
              #(listOfrowsByFeature[0].__getattribute__(opts.featName),
               #list(chrms),
               #opts.chrm)
    
    # --- Check strand info ---
    if opts.strand:
        strands = set([x.__getattribute__(opts.strand) for x in listOfrowsByFeature])
        
        # ensure all blocks on same strand
        if len(strands) != 1: 
            print >> sys.stderr, 'ERROR: the %s feature has more than one strand type: %s (%s column).  SKIPPING' % \
                  (listOfrowsByFeature[0].__getattribute__(opts.featName),
                   list(strands),
                   opts.strand)
            return None
            #raise Exception, \
                  #'ERROR: the %s feature has more than one strand type: %s (%s column)' % \
                  #(listOfrowsByFeature[0].__getattribute__(opts.featName),
                   #list(strands),
                   #opts.strand)
    
        # ensure the strand representation is recognized
        if not list(strands)[0] in strandReps: 
            raise Exception, \
                  "ERROR: the %s feature's strand representation is not recognized: %s (%s column)" % \
                  (listOfrowsByFeature[0].__getattribute__(opts.featName),
                   list(strands)[0],
                   opts.strand)
        
    # --- Check ThickEdge info ---
    # ensure all thkStarts agree
    if opts.thkStrt:
        thkStrts = set([x.__getattribute__(opts.thkStrt) for x in listOfrowsByFeature])
        if len(thkStrts) != 1: 
            raise Exception, \
                  'ERROR: the %s feature has more than one thickStart: %s (%s column)' % \
                  (listOfrowsByFeature[0].__getattribute__(opts.featName),
                   list(thkStrts),
                   opts.thkStrt)
    # ensure all thkEnds agree
    if opts.thkEnd:
        thkEnds = set([x.__getattribute__(opts.thkEnd) for x in listOfrowsByFeature])
        if len(thkEnds) != 1: 
            raise Exception, \
                  'ERROR: the %s feature has more than one thickEnd: %s (%s column)' % \
                  (listOfrowsByFeature[0].__getattribute__(opts.featName),
                   list(thkEnds),
                   opts.thkEnd)
    
    # +++ Sort based on blockStarts +++
    listOfrowsByFeature.sort(key=lambda x: x.__getattribute__(opts.blkChmStrt))
    
    # +++ Gather the required pieces +++    
    chrm      = listOfrowsByFeature[0].__getattribute__(opts.chrm)
    chrmStart = str(int(listOfrowsByFeature[0].__getattribute__(opts.blkChmStrt))-1)
    chrmEnd   = listOfrowsByFeature[-1].__getattribute__(opts.blkChmEnd)
    name      = listOfrowsByFeature[0].__getattribute__(opts.featName)
    score     = '0'
    rgb       = opts.rgb
    blkCount  = str(len(listOfrowsByFeature))
    blkSizes  = getBlockSizes(listOfrowsByFeature)
    blkStarts = getBlockStarts(listOfrowsByFeature,chrmStart)
    # --- Initialize pieces that depend on user options ---
    if opts.strand:
        strand = strandReps[listOfrowsByFeature[0].__getattribute__(opts.strand)]
    else:
        strand = '+' 
    if opts.thkStrt:
        thkStart  = listOfrowsByFeature[0].__getattribute__(opts.thkStrt)
    else:
        thkStart  = chrmStart
    if opts.thkEnd:
        thkEnd    = listOfrowsByFeature[0].__getattribute__(opts.thkEnd)
    else:
        thkEnd    = chrmEnd
    
    print '%s' % ('\t'.join([chrm,    
                             chrmStart,
                             chrmEnd,
                             name,
                             score,
                             strand,
                             thkStart,
                             thkEnd,
                             rgb,
                             blkCount,
                             blkSizes,
                             blkStarts]))


if __name__ == "__main__":
    
    
    #+++++++++++ File Parseing Etc +++++++++++    
    usage = """python %prog inFile [options]"""
    parser = optparse.OptionParser(usage)
    parser.add_option('-t',dest="track_name",type="str", default="untitled", 
                      help="""Space-less string for the display name for this track. (default=%default)""")
    parser.add_option('-d',dest="description",type="str", default="no description", 
                      help="""Quoted string: Exp: "Clone Paired Reads". (default=%default)""")
    parser.add_option('--featName',dest="featName",type="str", default=None, 
                      help="""REQUIRED - Exact Title of Column holding the names of the features. Exp: transcript_name (default=%default)""")
    parser.add_option('--blkChmStrt',dest="blkChmStrt",type="str", default=None, 
                      help="""REQUIRED - Exact Title of Column holding the chrom start positions of the blocks. Exp: seq_region_start (default=%default)""")
    parser.add_option('--blkChmEnd',dest="blkChmEnd",type="str", default=None, 
                      help="""REQUIRED - Exact Title of Column holding the chrom end positions of the blocks. Exp: seq_region_end (default=%default)""")
    parser.add_option('--chrm',dest="chrm",type="str", default=None, 
                      help="""REQUIRED - Exact Title of Column holding the name of the chrom on which each block resides. Exp: name -or- Chromosome/plasmid (default=%default)""")
    parser.add_option('--strand',dest="strand",type="str", default=None, 
                      help="""Exact Title of Column holding the strand indicator. Exp: strand. If None: all will be assigned '+' strand. (default=%default)""")
    parser.add_option('--thkStrt',dest="thkStrt",type="str", default=None, 
                      help="""Exact Title of Column holding the position of the thickStart. Exp: start_of_coding (default=%default)""")
    parser.add_option('--thkEnd',dest="thkEnd",type="str", default=None, 
                      help="""Exact Title of Column holding the position of the thickEnd. Exp: end_of_coding (default=%default)""")
    parser.add_option('--rgb',dest="rgb",type="str", default="0,0,0", 
                      help="""Space-less string to assign the color for this track. Exp: 0,0,0=black; 255,0,0=red; 0,255,0=green; 0,0,255=blue (default=%default)""")

    
    (opts, args) = parser.parse_args()
    
    if len(args) != 1:
        parser.print_help()
        exit()
    
    
    featBlkTable   = tableFile2namedTuple(args[0])
    rowsByFeature = groupFeatBlks(featBlkTable,opts)
    
    print """track name=%s description="%s" useScore=0""" % (opts.track_name, opts.description)
    
    for feature in rowsByFeature:
        printBEDline(rowsByFeature[feature],opts)
    