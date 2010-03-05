from time import time
from time import ctime
import optparse
import sys
import os
import scipy.stats as stats

try:
    import gusPyCode
except:
    sys.path.extend(['/Users/biggus/lib/python/site-packages/',
                     '/Library/Python/2.5/site-packages/',
                     '/Library/Frameworks/Python.framework/Versions/2.5/lib/python2.5/site-packages/'])
from gusPyCode.defs.JamesDefs import Bag
from gusPyCode.defs.bioDefs import ParseFastQ,ParseBowtieBed,ParseSolexaSorted,ParseBowtieMap


 

#+++++++++ Definitions ++++++++++ 

supportedFileTypes = {"bowtie_bed": ParseBowtieBed,
                      "bowtie_map": ParseBowtieMap,
                      "fastq"     : ParseFastQ,
                      "sorted"    : ParseSolexaSorted,}
def setParser(filePath,fileType):
    """Return the correct parser initialized and ready for use."""
    assert fileType in supportedFileTypes,\
           "** ERROR: %s not supported file type: %s **" % (fileType,supportedFileTypes.keys())
    return supportedFileTypes[fileType](filePath)

    
def addRead(readEntry,readDict,whichFile):
    """Increment count for read in correct file count. If not in dict,
    initialize entry and increment correct file's count."""
    
    if readEntry in readDict:
        readDict[readEntry][whichFile] += 1
    else:
        readDict[readEntry] = [0,0]
        readDict[readEntry][whichFile] += 1


if __name__ == '__main__':
    print '\n\n\n'
    
    #+++++++++ Parse Command Line ++++++++++
        
    usage = """python %prog -o outName readsFile1 readsFile2"""
    parser = optparse.OptionParser(usage)
    parser.add_option('-o',dest="out_name",type="string",default=False,
                      help="""<required> Name to give result files.""")
    parser.add_option('-t',dest="compare_type",type="string",default='readSeq',
                      help="""Type of data to compare [readSeq,readCoords] (default=%default)""")
    parser.add_option('-f',dest="file_type",type="string",default='fastq',
                      help="""File format of readsFiles %s (default=%s)""" % (supportedFileTypes.keys(),"%default"))
    parser.add_option('--split',dest="split",action="store_true",default=False,
                      help="""Attempt to split one file into two data sets based on the machine+lane information. (default=%default)""")
    parser.add_option('--names',dest="names",type="string",default=False,
                      help="""If using '--split' you _MUST_ give comma separated names to attach to datasets. Exp: --names first_dataSet,other_dataSet  (default=%default)""")


    
    
    (opts, args) = parser.parse_args()
    
    if len(sys.argv) == 1:
        parser.print_help()
        exit(1)
    if (len(args) != 2) and (opts.split == False):
        parser.print_help()
        print "\n\n** ERROR: Please supply exactly two readsFiles or use the '--split' option. **"
        exit(1)
    if (opts.file_type not in ['bowtie_bed', 'bowtie_map']) and (opts.split == True):
        parser.print_help()
        print "\n\n** ERROR: The '--split' option is ONLY allowed with file_type(bowtie_bed). **"
        exit(1)
    if (opts.split == True) and (opts.names == False):
        parser.print_help()
        print "\n\n** ERROR: When using '--split' you _MUST_ use '--names'. **"
        exit(1)
    if not opts.out_name:
        parser.print_help()
        print "\n\n** ERROR: You must supply an out name. **"
    if opts.compare_type not in ['readSeq','readCoords']:
        parser.print_help()
        print "\n\n** ERROR: compare_type can only be %s **" % (['readSeq','readCoords'])
    
    readCoordIncompatibles = ['fastq',]
    if (opts.compare_type  == 'readCoords') and (opts.file_type in readCoordIncompatibles):
        parser.print_help()
        print "\n\n** ERROR: compare_type(%s) is not compatible with file_types(%s)**" % (opts.compare_type,readCoordIncompatibles)
        

    # ++++++++++ Open Files And Get Variables Set up ++++++++
    readDict  = {} ## to collect read counts
    whichFile = None
    ## ++++ if only one file: ++++
    if opts.split:
        oneFile = args[0]
        # -- set up correct parser --
        oneParser = setParser(oneFile,opts.file_type)
        if opts.compare_type == 'readCoords':
            parseEntryType = oneParser._parseCoords
        elif opts.compare_type == 'readSeq':
            parseEntryType = oneParser._parseReadSeq
        
            
        # ++++++++++ Execute the counting ++++++++
        print 'Counting single file...'
        t1_0 = time()
        whichRep = 0
        readLine = oneParser.getNext()
        while 1:
            lastLine = readLine[:]
            readEntry = parseEntryType(readLine)
            if readEntry:
                addRead(readEntry,readDict,whichRep)
            else: break
            # --- if moved to the next dataset, increment whichRep ---
            readLine = oneParser.getNext()
            if readLine == None:
                break
            else:
                if readLine[0].split(':')[:3] != lastLine[0].split(':')[:3]:
                    whichRep += 1
                    assert whichRep <= 1, \
                           """** ERROR: It looks like your bowtie.map file has more than one dataset
                           based on the machine+lane criteria I use. **"""
                        
                
        t2_0 = time()
        
    ## ++++ if two files: ++++    
    else:
        readsFile0   = args[0]
        readsFile1   = args[1]
        # -- set up correct parsers --
        parser0 = setParser(readsFile0,opts.file_type)
        parser1 = setParser(readsFile1,opts.file_type)
        if opts.compare_type == 'readCoords':
            parser0_getNext = parser0.getNextReadCoords
            parser1_getNext = parser1.getNextReadCoords
        elif opts.compare_type == 'readSeq':
            parser0_getNext = parser0.getNextReadSeq
            parser1_getNext = parser1.getNextReadSeq
    
        
        
        
        
        # ++++++++++ Execute the counting ++++++++
        # --- readsFile0 ---
        print 'Counting the first file...'
        t1_0 = time()
        whichRep = 0
        while 1:
            readEntry = parser0_getNext() 
            if readEntry:
                addRead(readEntry,readDict,whichRep)
            else: break
        t2_0 = time()
        print 'Counting took %s min.' % ((t2_0-t1_0)/60)
        # --- readsFile1 ---
        print 'Counting the second file...'
        t1_1 = time()
        whichRep = 1
        while 1:
            readEntry = parser1_getNext()
            if readEntry:
                addRead(readEntry,readDict,whichRep)
            else: break
        t2_1 = time()
        print 'Counting took %s min.' % ((t2_1-t1_1)/60)
        
    # ++++++++++ Calculate The Correlations ++++++++

    assert whichRep == 1, \
           """** ERROR: It looks like your bowtie.map file has ONLY one dataset
           based on the machine+lane criteria I use. **"""

    
    print 'Getting vectors... %s' % (ctime())
    vector0 = tuple([x[0] for x in readDict.values()])
    vector1 = tuple([x[1] for x in readDict.values()])
    print ctime()
    
    print 'Deleting readDict...'
    t1_del = time()
    del(readDict)
    t2_del = time()
    print '%s min.' % ((t2_del-t1_del)/60)
    
    print 'Calculating Pearson... %s' % (ctime())
    pearson  = stats.pearsonr(vector0,vector1)
    print pearson

    
    # ++++++++++ write labeled vectors to file ++++++++ 
    print 'Writing labeled vectors to file....'
    oFile = open(opts.out_name,'w')
    if opts.names:
        oFile.write('%s\t%s\n' % (opts.names.split(',')[0],opts.names.split(',')[1])) # label columns
    else:
        oFile.write('%s\t%s\n' % (readsFile0.split('/')[-1],readsFile1.split('/')[-1])) # label columns
    for i in range(len(vector0)):
        oFile.write('%s\t%s\n' % (vector0[i],vector1[i]))
    oFile.flush()
    oFile.close()