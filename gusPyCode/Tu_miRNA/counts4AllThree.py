from gusPyCode.defs.JamesDefs import odd_or_even
import re

noteToSelf = \
           """
           So far I have set this up to count miRNAs per gene.  I need to count genes per miRNA.
           may be able to salvage this work using the zip(*list) trick to transpose the data matrices.
           """


def processFiles(pathList):
    """
    Load, group pos1 and pos2 data, and return dict of dictionaries of hits for match/ctrl
    source files using genome prefix as keys for master dict.
    """
    lettersRegEx = re.compile('^\D+', re.IGNORECASE)
    files = []
    masterList = []
    
    for i in range(len(pathList)):
        data = map(lambda l: l.strip().split('\t'), open(pathList[i],'rU'))
        if data[0][0].startswith('#'): data.pop(0) # get rid of headers
        tmpDict = {}
        for i in range(len(data)):
            tmpDict[data[i][0]] = combinePos1and2(data[i][1:])
        
        
        # Add new hitDict to the master list
        masterList.append(tmpDict)
    return masterList

def combinePos1and2(listOfhits):
    pos1s = []
    pos2s = []
    for i in range(len(listOfhits)):
        if odd_or_even(i) == 'even':  # remember 0 == even
            pos1s.append(eval(listOfhits[i]))
        else:
            pos2s.append(eval(listOfhits[i]))
    assert len(pos1s) == len(pos2s), \
           'Error: len(pos1s) _MUST_ == len(pos2s)'
    
    # zip trick them back together and return the list of lists
    rList = zip(*[pos1s,pos2s])
    return rList

def filter4genes(processedData):
    pass


threeFile  = '/Users/biggus/Documents/James/Data/OrthologDefs/nrOrthos/!nrOrthoDefs/3wayOrthos.combineOrthosJamesDefs.out.filtered.txt'

# MUST insure that matchFile indeces match the ctrlFile indeces!!!
matchFiles = ['/Users/biggus/Documents/James/Data/Tu_miRNA/SeedCountOutPut/counts/miRBaseMatureSeedsOn_Aa_500afterCoding.match.txt',
              '/Users/biggus/Documents/James/Data/Tu_miRNA/SeedCountOutPut/counts/miRBaseMatureSeedsOn_Ag_500afterCoding.match.txt',
              '/Users/biggus/Documents/James/Data/Tu_miRNA/SeedCountOutPut/counts/miRBaseMatureSeedsOn_Cq_500afterCoding.match.txt']

ctrlFiles  = ['/Users/biggus/Documents/James/Data/Tu_miRNA/SeedCountOutPut/counts/miRBaseMatureSeedsOn_Aa_500afterCoding.ctrl.txt',
              '/Users/biggus/Documents/James/Data/Tu_miRNA/SeedCountOutPut/counts/miRBaseMatureSeedsOn_Ag_500afterCoding.ctrl.txt',
              '/Users/biggus/Documents/James/Data/Tu_miRNA/SeedCountOutPut/counts/miRBaseMatureSeedsOn_Cq_500afterCoding.ctrl.txt']



defsAll3 = map(lambda l: l.strip(), open(threeFile, 'rU'))
defsAll3 = ','.join(defsAll3)
defsAll3 = eval(defsAll3)

matchData = processFiles(matchFiles)
ctrlData  = processFiles(ctrlData)

matchesInOne   = [0,0]
matchesInTwo   = [0,0]
matchesInThree = [0,0]

ctrlsInOne   = [0,0]
ctrlsInTwo   = [0,0]
ctrlsInThree = [0,0]


ltrRE = re.compile('^\D+', re.IGNORECASE)

