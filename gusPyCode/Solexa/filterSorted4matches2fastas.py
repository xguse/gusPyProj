from time import time
import optparse
import sys
import motility
from TAMO.seq import Fasta

def fastas2Dict(listOfPaths):
    """ """
    rFastaDict = {}
    for path in listOfPaths:
        tempDict = Fasta.file2dict(path)
        rFastaDict.update(tempDict)
    return rFastaDict

def iterReadsLine(readsFileLine,newReadsFile,fastaDict,hitDict,fileType):
    """ """
    fileTypes = {"sorted":8,
                 "bowtie":4}
    if not fileType == "bowtie":
        print "** ERROR: At this time only 'bowtie' fileType is supported **"
        exit(1)  
    if not fileType in fileTypes:
        parser.print_help()
        print "** ERROR: %s is not a valid fileType: %s **" % (opts.file_type, fileTypes.keys())
        exit(1)       
    line = readsFileLine
    read = line.split('\t')[fileTypes[opts.file_type]]
    searchResult = searchFastas(read,fastaDict)
    if len(searchResult) == 0:
        writeOut(line.split('\t'),newReadsFile,inType=opts.file_type,outType=opts.out_type)
    else:
        for result in list(searchResult):
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
    

def searchFastas(read,fastaDict):
    """ """
    tmpSeqNames = set()
    for seqName in fastaDict:
        matches = motility.find_iupac(fastaDict[seqName], read, mismatches=2)
        if not len(matches) == 0:
            tmpSeqNames.add(seqName)
    return tmpSeqNames
    
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
    parser.add_option('-f',dest="fasta_files",type="string",default=False,
                      help="""<required> Comma separated list (no spaces) of fasta files to use as search material.""")
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
    if not opts.fasta_files:
        parser.print_help()
        print "\n\n** ERROR: Please supply at least one fasta file. **"
    if not opts.run_name:
        parser.print_help()
        print "\n\n** ERROR: You must supply a run name. **"
        
    opts.fasta_files = opts.fasta_files.split(',')

    # ++++++++++ Open Files And Get Variables Set up ++++++++
    readsFile    = open(args[0],'rU')
    newReadsFile = open("%s_filtered.txt" % (opts.run_name),'w')
    fastaDict    = fastas2Dict(opts.fasta_files)
    hitDict      = {} ## to collect reads that match the fasta
    
    # -- init hitDict --
    for name in fastaDict:
        hitDict[name] = 0
        
    # ++++++++++ Its Bidness Time! ++++++++
    t1 = time() ## Start timer
    while 1:
        readsFileLine = readsFile.readline().strip('\n')
        if readsFileLine:
            iterReadsLine(readsFileLine,newReadsFile,fastaDict,hitDict,opts.file_type)
        else: break
    
    # ++++++++++ Write OutTable ++++++++
    outTable = open("%s_hitTable.txt" % (opts.run_name), 'w')
    for seqName in sorted(hitDict):
        if hitDict[seqName] > 0:
            outTable.write('%s\t%s\n' % (seqName, hitDict[seqName]))
    
    t2 = time() ## End timer    
    print 'This run took %s min.' % ((t2-float(t1))/60)
    
    # ++++++++++ Close Up Shop ++++++++
    readsFile.close()
    newReadsFile.close()
    outTable.close()