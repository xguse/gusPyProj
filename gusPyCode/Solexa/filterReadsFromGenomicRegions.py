from time import time
import optparse
import sys
from gusPyCode.defs.JamesDefs import detect_1D_overlap
from gusPyCode.defs.JamesDefs import Bag

def any2noChange(lineList):
    """Returns the lineList as is."""
    return lineList
    
def sorted2DEGseq(lineList):
    """Code to convert sorted to DEGseq."""
    rList = [str(lineList[11]),
             str(lineList[12]),
             str(int(lineList[12])+int(lineList[14])-1),
             '|'.join(lineList[:8])+'_'+lineList[8],
             str(0),
             lineList[13]]
    if   lineList[-1] == 'R':lineList[-1]='-'
    elif lineList[-1] == 'F':lineList[-1]='+'
    
    return rList
    
def bowtie2DEGseq(lineList):
    """Code to convert bowtie to DEGseq."""
    return [lineList[2],
            lineList[3],
            str(len(lineList[4])+int(lineList[3])-1),
            lineList[0]+'_'+lineList[4],
            lineList[6],
            lineList[1]]

def getCoords(filePath):
    """Get and parse list of "restricted" genomic alignment coords."""
    return map(lambda l: l.strip('\n').split('\t'), open(filePath,'rU'))

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
            
    

def iterReadsLine(readsFileLine,newReadsFile,coordsDict,hitDict,fileType):
    """Call checkCoords() for overlap in coords and write to filtered file."""
    fileTypes = {"sorted":Bag({'seqName':'line[11]',
                               'start':'int(line[12])',
                               'end':'int(line[12])+int(line[14])-1'}),
                 "bowtie":Bag({'seqName':'line[2]',
                               'start':'int(line[3])',
                               'end':'int(line[3])+len(line[4])-1'})}
 
    if not fileType in fileTypes:
        parser.print_help()
        print "** ERROR: %s is not a valid fileType: %s **" % (opts.file_type, fileTypes.keys())
        exit(1)       
    line = readsFileLine.split('\t')
    fileStats = fileTypes[opts.file_type]
    readCoords = [eval(fileStats.seqName), eval(fileStats.start), eval(fileStats.end)]
    searchResult = checkCoords(readCoords,coordsDict)
    if not searchResult:
        writeOut(line,newReadsFile,inType=opts.file_type,outType=opts.out_type)
    else:
        hitDict[searchResult] += 1
            
    

def checkCoords(readCoords,coordsDict):
    """Check for overlap in coords"""
    
    hitDictKey = False
    # -- get right contig data if exists --
    if coordsDict.has_key(readCoords[0]):
        cntgNixCoords = coordsDict[readCoords[0]]
    else:
        return hitDictKey
    # -- cycle through nixCoords for overlap --
    for i in range(len(cntgNixCoords)):
        if detect_1D_overlap(readCoords[1:3],cntgNixCoords[i][0:2]):
            # -- If hit, set hitDictKey and break --
            hitDictKey = tuple([readCoords[0]]+cntgNixCoords[i]) # re-constructs a key in hitDict
            break
    

    return hitDictKey
    
def writeOut(line,newReadsFile,inType="bowtie",outType="DEGseq"):
    """Write line to new file in correct file format."""
    cvrtTo     = ("no_change", 
                  "DEGseq")
    
    cvrtFrom   = ("sorted", 
                  "bowtie")
    
    cvrtnTable = {"bowtie":{"DEGseq"    : bowtie2DEGseq,
                            "no_change" : any2noChange,},
                  "sorted":{"DEGseq"    : sorted2DEGseq,
                            "no_change" : any2noChange,},
                  }
    
    if not inType in cvrtFrom:
        parser.print_help()
        print "** ERROR: %s is not a valid inType: %s **" % (inType, cvrtFrom)
        exit(1)  
    if not outType in cvrtTo:
        parser.print_help()
        print "** ERROR: %s is not a valid outType: %s **" % (outType, cvrtTo)
        exit(1)
        
    toWrite = cvrtnTable[inType][outType](line)
    newReadsFile.write('\t'.join(toWrite+['\n']))
    


if __name__ == '__main__':
    print '\n\n\n'
    
    #+++++++++ Parse Command Line ++++++++++
        
    usage = """python %prog readsFile -g genomicCoordsFile -r runName [options]"""
    parser = optparse.OptionParser(usage)
    parser.add_option('-g',dest="genomic_coords",type="string",default=False,
                      help="""<required> A filter file containing lines of form:\nseqName<tab>start<tab>stop""")
    parser.add_option('-r',dest="run_name",type="string",default=False,
                      help="""<required> Base name to give result files.""")
    parser.add_option('-t',dest="file_type",type="string",default="bowtie",
                      help="""Which type of file is readsFile? ["sorted","bowtie"] (default=%default)""")
    parser.add_option('-o',dest="out_type",type="string",default="no_change",
                      help="""Write new file in format ["DEGseq","no_change"] (default=%default)""")
    parser.add_option('-c',dest="count",action="store_true",default=False,
                      help="""Only print the count of reads that matched the fastas. (default=%default)""")
    
    
    (opts, args) = parser.parse_args()
    
    if len(sys.argv) == 1:
        parser.print_help()
        exit(1)
    if len(args) == 0:
        parser.print_help()
        print "\n\n** ERROR: Please supply readsFile. **"
        exit(1)
    if not opts.genomic_coords:
        parser.print_help()
        print "\n\n** ERROR: Please supply at least one fasta file. **"
    if not opts.run_name:
        parser.print_help()
        print "\n\n** ERROR: You must supply a run name. **"
        


    # ++++++++++ Open Files And Get Variables Set up ++++++++
    readsFile    = open(args[0],'rU')
    newReadsFile = open("%s_filtered.txt" % (opts.run_name),'w')
    coordsList   = getCoords(opts.genomic_coords)
    coordsDict   = buildCoordsDict(coordsList)
    hitDict      = {} ## to collect reads that match the fasta
    
    # -- init hitDict --
    for i in range(len(coordsList)):
        hitDict[tuple(coordsList[i])] = 0
        
    # ++++++++++ Its Bidness Time! ++++++++
    t1 = time() ## Start timer
    while 1:
        readsFileLine = readsFile.readline().strip('\n')
        if readsFileLine:
            iterReadsLine(readsFileLine,newReadsFile,coordsDict,hitDict,opts.file_type)
        else: break
    
    # ++++++++++ Write OutTable ++++++++
    outTable = open("%s_hitTable.txt" % (opts.run_name), 'w')
    for key in sorted(hitDict):
        if hitDict[key] > 0:
            outTable.write('%s\t%s\t%s\t%s\t%s\n' % (key[0],
                                                     key[1],
                                                     key[2],
                                                     key[3],
                                                     hitDict[key]))
    
    t2 = time() ## End timer    
    print 'This run took %s sec.' % (t2-float(t1))
    
    # ++++++++++ Close Up Shop ++++++++
    readsFile.close()
    newReadsFile.close()
    outTable.close()