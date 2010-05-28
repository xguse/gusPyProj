import sys
import csv
import optparse
from gusPyCode.defs.JamesDefs import filter2Dlist
from gusPyCode.defs.mathDefs import mean


# ++++ Definitions ++++
def averageProbeSets(xpnData):
    """return xpnData after averageing probeset results for genes with
    multiple probesets."""
    newXpnData = []
    dataDict   = {}
    
    for line in xpnData:
        if line[1] in dataDict:
            dataDict[line[1]].append(line)
        else:
            dataDict[line[1]] = []
            dataDict[line[1]].append(line)
            
    for gene in dataDict:
        newDataLine = ["avg'd",gene]
        for i in range(2,len(dataDict[gene][0])):
            newDataLine.append(mean([float(x) for x in [y[i] for y in dataDict[gene]]]))
        
        newXpnData.append(tuple(newDataLine[:]))
    
    return newXpnData
    
def getExpressionData(path2File):
    """Load data from mircroarray."""
    maFile = open(path2File,'rU')
    reader = csv.reader(maFile)
    
    xpnData = []
    for line in reader:
        xpnData.append(line)
    headers = xpnData.pop(0)
    
    return (headers,xpnData)

def getMicroRNAData(path2File):
    """Load miRNA target prediction file into a list of rows and fields"""
    mirFile = open(path2File, 'rU')
    
    return map(lambda l: l.strip('\n').split(' : '), open(path2File, 'rU'))
    
    
def cycleTimePoints(xpnData,tpCol,headers,mirData,cnvClass,seedType):
    """Given the timepoint, sort it and cycle through the relevant miR/target sets
    to writeout the 0s and 1s for genes."""
    
    # ++++ Sort xpnData ++++
    xpnData.sort(key=lambda x: float(x[tpCol]))
    
    # ++++ Filter Desired mirTarget Sets ++++
    if seedType.startswith('all'): 
        trgtGroups = filter2Dlist(mirData,1,'allPassedSeedsFor_%s' % (cnvClass.split('_')[-1]))  # filter orthoType
    else:
        trgtGroups = filter2Dlist(mirData,2,cnvClass)     # filter orthoType
        trgtGroups = filter2Dlist(trgtGroups,1,seedType)  # filter seedType
    
    for i in range(len(trgtGroups)):
        writeOut0and1s(xpnData, trgtGroups[i],headers,tpCol,cnvClass,seedType)
    

def writeOut0and1s(xpnData,mirTrgtList,headers,tpCol,cnvClass,seedType):
    """For a miR-seedType: write the 0sAnd1s to a
    file named -> miR.seedType.orthoType.timePoint.0s1s.txt"""
    
    if seedType.startswith('all'):
        outName = '%s.%s.%s.0s1s.txt' % (mirTrgtList[0],
                                         mirTrgtList[1],
                                         headers[tpCol])
    else:
        outName = '%s.%s.%s.%s.0s1s.txt' % (mirTrgtList[0],
                                            mirTrgtList[1],
                                            mirTrgtList[2],
                                            headers[tpCol])
    oFile = open(outName,'w')
    
    mirTrgtGenes = getAGAP(mirTrgtList[-1])
    rankedGenes = tuple([x[1] for x in xpnData])
    for gene in rankedGenes:
        if gene in mirTrgtGenes:
            oFile.write('1\n')
        else:
            oFile.write('0\n')
            
    oFile.flush()
    oFile.close()

def getAGAP(strOfTupReps):
    rSet = set()
    tups = eval(strOfTupReps)
    
    for tup in tups:
        for i in tup:
            if i.startswith('AGAP'):
                rSet.add(i)
    return rSet
    

seedTypes = ['m3_to_m8',
             'm2_to_m7',
             'A1_to_m7',
             'm2_to_m8',
             'A1_to_m8',
             'all_passed']
    
if __name__ == '__main__':
    
    print '\n'
    
    # +++++ Parse Command line +++++
    usage = """python %prog [options] microArrayCSV mirTargetFile"""
    parser = optparse.OptionParser(usage)

    parser.add_option('-s',dest="seed_type",type="string",default='all_passed',
                      help="""Which SeedType to use. (default=%default)""")
    parser.add_option('-o',dest="ortho_type",type="int",default=3,
                      help="""Which orthoType to use [2,3]. (default=%default)""")
    parser.add_option('-t',dest="time_point",type="string",default='all',
                      help="""Which timePoint column to use [3,4,5,6,7,all]. (default=%default)""")

    (opts, args) = parser.parse_args()
    
    if len(sys.argv) == 1:
        parser.print_help()
        exit(1)
    if len(args) < 2:
        parser.print_help()
        print "\n\n** ERROR: Please supply both microArrayCSV and mirTargetFile. **"
        exit(1)
    if opts.ortho_type not in [2,3]:
        parser.print_help()
        print "\n\n** ERROR: -o may only be [2,3]. **" 
        exit(1)
    if opts.seed_type not in seedTypes:
        parser.print_help()
        print "\n\n** ERROR: -s may only be %s. **" % (seedTypes)
        exit(1)
    if opts.time_point not in ['3','4','5','6','7','all']:
        parser.print_help()
        print "\n\n** ERROR: -t may only be ['3','4','5','6','7','all']. **"
        exit(1)
        
    
    print 'Reading in data files...'
    headers,xpnData = getExpressionData(args[0])
    xpnData         = averageProbeSets(xpnData)
    mirData         = getMicroRNAData(args[1])
    

    # +++++ Set up the Cycle Routine +++++
    # --- what time points are we doing ---
    if opts.time_point != 'all':
        tpCols = [int(opts.time_point)]
    else:
        tpCols = range(3,8)
    # --- intialize other vars ---    
    cnvClass = 'orthoType_%s' % (opts.ortho_type)
    seedType = opts.seed_type
    
    # +++++ run Cycle +++++
    for tpCol in tpCols:
        cycleTimePoints(xpnData,tpCol,headers,mirData,cnvClass,seedType)
