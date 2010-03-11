from time import time
import optparse
import sys
from gusPyCode.defs.JamesDefs import detect_1D_overlap
from gusPyCode.defs.JamesDefs import Bag
from gusPyCode.defs.mosqData.AaCntgCvrt import supContigConvert


def getCoords(filePath):
    """Get and parse list of "restricted" genomic alignment coords."""
    coordsFile = open(filePath,'rU')
    coordsList = []
    for line in coordsFile:
        tmp = line.strip('\n').split('\t')
        coordsList.append([supContigConvert[tmp[0]]]+tmp[1:])
        
    return coordsList
        

def buildCoordsDict(coordsList):
    """Kludge to build dict allowing for faster lookups while not breaking what "needs" coordsDict.
    keys = contigName; values = listOfCoordsOnContig"""
    
    coordsDict = {}
    for i in range(len(coordsList)):
        if coordsDict.has_key(coordsList[i][0]):
            coordsDict[coordsList[i][0]].append(coordsList[i][1:])
        else:
            coordsDict[coordsList[i][0]] = []
            coordsDict[coordsList[i][0]].append(coordsList[i][1:])
            
    return coordsDict
            
    

def iterSNPsLine(snpsFileLine,coordsDict,hitDict,orphans):
    """Call checkCoords() for overlap in coords and write to filtered file."""
    
    if int(snpsFileLine[4]) <= opts.coverage:
        return
    
    readCoords = [snpsFileLine[0], snpsFileLine[1], snpsFileLine[1]]
    searchResult = checkCoords(readCoords,coordsDict)
    if not searchResult:
        orphans += 1
    else:
        hitDict[searchResult].append(snpsFileLine)
            
    

def checkCoords(snpCoords,coordsDict):
    """Check for overlap in coords"""
    
    hitDictKey = False
    # -- get right contig data if exists --
    if coordsDict.has_key(snpCoords[0]):
        genCoords = coordsDict[snpCoords[0]]
    else:
        return hitDictKey
    # -- cycle through genCoords for overlap --
    for i in range(len(genCoords)):
        if detect_1D_overlap(snpCoords[1:3],genCoords[i][0:2]):
            # -- If hit, set hitDictKey and break --
            hitDictKey = tuple([snpCoords[0]]+genCoords[i]) # re-constructs a key in hitDict
            break
    

    return hitDictKey
    
    


if __name__ == '__main__':
    print '\n\n\n'
    
    #+++++++++ Parse Command Line ++++++++++
        
    usage = """python %prog snpsFile -g genomicCoordsFile -r runName [options]"""
    parser = optparse.OptionParser(usage)
    parser.add_option('-g',dest="genomic_coords",type="string",default=False,
                      help="""<required> A filter file containing lines of form:\nseqFrag<tab>start<tab>stop<tab>name""")
    parser.add_option('-r',dest="run_name",type="string",default=False,
                      help="""<required> Base name to give result files.""")
    parser.add_option('--cov',dest="coverage",type="int",default=10,
                      help="""Required read coverage at SNP site. (default=%default)""")
    
    (opts, args) = parser.parse_args()
    
    if len(sys.argv) == 1:
        parser.print_help()
        exit(1)
    if len(args) == 0:
        parser.print_help()
        print "\n\n** ERROR: Please supply snpsFile. **"
        exit(1)
    if not opts.genomic_coords:
        parser.print_help()
        print "\n\n** ERROR: Please supply genomicCoordsFile. **"
    if not opts.run_name:
        parser.print_help()
        print "\n\n** ERROR: You must supply a run name. **"
        


    # ++++++++++ Open Files And Get Variables Set up ++++++++
    snpsFile   = open(args[0],'rU')
    coordsList = getCoords(opts.genomic_coords)
    coordsDict = buildCoordsDict(coordsList)
    hitDict    = {} ## to collect reads that match the GenCoords
    orphans    = 0  ## to count the number of SNPs not mapping to known exons

    
    # -- init hitDict --
    for i in range(len(coordsList)):
        hitDict[tuple(coordsList[i])] = []
        
    # ++++++++++ Its Bidness Time! ++++++++
    t1 = time() ## Start timer
    snpLine = 0
    while 1:
        snpsFileLine = snpsFile.readline()
        if snpsFileLine:
            snpLine += 1
            snpsFileLine = snpsFileLine.strip('\n').split('\t')
            try:
                int(snpsFileLine[1])
            except:
                continue
            iterSNPsLine(snpsFileLine,coordsDict,hitDict,orphans)
            if snpLine%5000 == 0:
                print "snpLine: %s" % (snpLine)
        else: 
            break
    
    # ++++++++++ Write OutTable ++++++++
    outTable = open("%s.snpsPerBp.cov%s.txt" % (opts.run_name,opts.coverage), 'w')
    outTable.write('# orphan snps: %s\n' % (orphans))
    for key in sorted(hitDict):
        if hitDict[key] > 0:
            contig       = key[0]
            exonName     = key[3]
            start        = key[1]
            stop         = key[2]
            snps         = len(hitDict[key])
            snpsPerBp    = float(snps)/(int(stop)-int(start)+1)
            snpFileLines = str(hitDict[key])
            outTable.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (contig,
                                                             exonName,
                                                             start,
                                                             stop,
                                                             snps,
                                                             snpsPerBp,
                                                             snpFileLines))
    
    t2 = time() ## End timer    
    print 'This run took %s sec.' % (t2-float(t1))
    
    
    # ++++++++++ Close Up Shop ++++++++
    snpsFile.close()
    outTable.close()