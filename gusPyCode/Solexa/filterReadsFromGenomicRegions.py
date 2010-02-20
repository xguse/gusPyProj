from time import time
import optparse
import sys
from TAMO.seq import Fasta
from gusPyCode.defs.JamesDefs import Bag

def getCoords(filePath):
    """ """
    return map(lambda l: l.strip('\n').split('\t'), open(filePath,'rU'))

def iterReadsLine(readsFileLine,newReadsFile,coordsList,hitDict,fileType):
    """ """
    fileTypes = {"sorted":{},
                 "bowtie":Bag({'seqName':'line[2]',
                               'start':'int(line[3])',
                               'end':'int(line[3])+len(line[4])-1'})}
    if not fileType == "bowtie":
        print "** ERROR: At this time only 'bowtie' fileType is supported **"
        exit(1)  
    if not fileType in fileTypes:
        parser.print_help()
        print "** ERROR: %s is not a valid fileType: %s **" % (opts.file_type, fileTypes.keys())
        exit(1)       
    line = readsFileLine.split('\t')
    fileStats = fileTypes[opts.file_type]
    readCoords = [eval(fileStats.seqName), eval(fileStats.start), eval(fileStats.end)]
    searchResult = checkCoords(readCoords,coordsList)
    if len(searchResult) == 0:
        writeOut(line,newReadsFile,inType=opts.file_type,outType=opts.out_type)
    else:
        for result in searchResult:
            hitDict[result] += 1
            
#def iterBowtie(readsFileLine,newReadsFile,fastaDict,hitDict):
    #""" """
    #line = readsFileLine
    #read = line.split('\t')[4]
    #searchResult = searchFastas(read,newReadsFile,fastaDict)
    #if len(searchResult) == 0:
        #newReadsFile.write(line)
    #else:
        #for result in searchResult:
            #hitDict[result] += 1
    

def checkCoords(readCoords,coordsList):
    """ """
    indexHits = []
    for i in range(len(coordsList)):
        # -- on right seqName? --
        if not readCoords[0] == coordsList[i][0]:
            continue
        # -- coords overlap with nixCoords? --
        lEdge = (int(readCoords[1])-int(coordsList[i][1]) > 0) and (int(readCoords[1])-int(coordsList[i][2]) < 0)
        rEdge = (int(readCoords[2])-int(coordsList[i][1]) > 0) and (int(readCoords[2])-int(coordsList[i][2]) < 0)
        # -- did we get a hit? --
        if lEdge or rEdge:
            indexHits.append(i)
            break
    return indexHits
    
def writeOut(line,newReadsFile,inType="bowtie",outType="DEGseq"):
    """ """
    cvrtTo     = ("no_change", 
                  "DEGseq")
    
    cvrtFrom   = ("sorted", 
                  "bowtie")
    
    cvrtnTable = {"bowtie":{"DEGseq":'''[line[2],line[3],str(len(line[4])+int(line[3])-1),line[0]+'_'+line[4],line[6],line[1]]''',
                            "no_change":'''line''',},
                  "sorted":{"DEGseq":'''exit('\n\n** ERROR: this conversion has not been defined yet. **\n** See def writeOut. **')''',
                            "no_change":'''line''',},
                  }
    
    if not inType in cvrtFrom:
        parser.print_help()
        print "** ERROR: %s is not a valid inType: %s **" % (inType, cvrtFrom)
        exit(1)  
    if not outType in cvrtTo:
        parser.print_help()
        print "** ERROR: %s is not a valid outType: %s **" % (outType, cvrtTo)
        exit(1)
        
    toWrite = eval(cvrtnTable[inType][outType])
    newReadsFile.write('\t'.join(toWrite+['\n']))
    


if __name__ == '__main__':
    print '\n\n\n'
    
    #+++++++++ Parse Command Line ++++++++++
        
    usage = """python %prog readsFile -f fasta1,fasta2,..fastaN -r runName """
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
    hitDict      = {} ## to collect reads that match the fasta
    
    # -- init hitDict --
    for i in range(len(coordsList)):
        hitDict[i] = 0
        
    # ++++++++++ Its Bidness Time! ++++++++
    t1 = time() ## Start timer
    while 1:
        readsFileLine = readsFile.readline().strip('\n')
        if readsFileLine:
            iterReadsLine(readsFileLine,newReadsFile,coordsList,hitDict,opts.file_type)
        else: break
    
    # ++++++++++ Write OutTable ++++++++
    outTable = open("%s_hitTable.txt" % (opts.run_name), 'w')
    for key in sorted(hitDict):
        if hitDict[key] > 0:
            outTable.write('%s\t%s\t%s\t%s\t\n' % (coordsList[key][0],
                                                     coordsList[key][1],
                                                     coordsList[key][2],
                                                     hitDict[key]))
    
    t2 = time() ## End timer    
    print 'This run took %s sec.' % (t2-float(t1))
    
    # ++++++++++ Close Up Shop ++++++++
    readsFile.close()
    newReadsFile.close()
    outTable.close()