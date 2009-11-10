from gusPyCode.defs import JamesDefs
from time import time

def getLittleK(countingDict, cluster, keysNotFound):
    k = 0
    
    for rec in cluster:
        if countingDict.has_key(rec[1]):
            k += countingDict[rec[1]]
        else:
            keysNotFound.add(rec[1]) 
    
    return k
    

def getBigK(countingDict):
    K = 0
    for AGAP in countingDict:
        K += countingDict[AGAP]
        
    return K

def findMotifComboInAGAP(AGAP, matrixOfAGAPsvMotifs, motifCombo, countingDict, motif2indexDict):
    
    #  Select first combo member after copying WITH OUT reference to new list
    motifComboCopy = motifCombo[:len(motifCombo)]
    comboMember = motifComboCopy.pop(0)
    
    #  test for combo member
    if matrixOfAGAPsvMotifs[AGAP][motif2indexDict[comboMember]] == 1:
        #  -check to see if this is last member
        if motifComboCopy == []:
            #  --if last member, place count in countingDict for AGAP
            countingDict[AGAP] = 1
            return
            
        else:
            #  --if not last member, call next instance of findMotifComboInAGAP
            findMotifComboInAGAP(AGAP, matrixOfAGAPsvMotifs, motifComboCopy, countingDict, motif2indexDict)
    else:
        #  -if ever a search fails stop searching and return to main
        return

def generateMotifCombos(motifList, mostMotifsInSet=3):
    from probstat import Combination
    
    listOfSearchCombos = []
    comboLen = 2
    lenOfmotifList = len(motifList)
    tl1 = time()
    while comboLen <= lenOfmotifList and comboLen <= mostMotifsInSet:
        comboObj = Combination(motifList,comboLen)
        for item in comboObj:
            listOfSearchCombos.append(item)
        comboLen+=1
    
    len_listOfSearchCombos = len(listOfSearchCombos)
    tl2 = time()
    
    print "%s motif combos produced: %s min" % (str(len_listOfSearchCombos),str((tl2-tl1)/60))
    return listOfSearchCombos


def readInMatrixFromFile(pathToFile):
    from time import time
    
    """
    Returns a list with:
    [0] a dict with AGAPs as keys and a list of 0s and 1s as the value for each AGAP
    [1] a list of motifs derived from the first line of the file which should have same order as respective 1s and 0s
    """
    t1 = time()
    # open and read in the file  
    fileList = map(lambda line : line.strip(), open(pathToFile, 'rU').readlines())
    
    #  create motifList from first line *** dont need to exclude first \t bc it was ignored in line above
    motifList = fileList.pop(0).split('\t')
    
    
    #  create dictOfAGAPS
    
    dictOfAGAPS = {}
    for line in fileList:
##        tAGAP1 = time()
        lineList = line.split('\t')
        AGAPstr = lineList.pop(0)
        
        #  convert 1s and 0s back to ints
        for i in range(len(lineList)):
            lineList[i] = int(lineList[i])

        dictOfAGAPS[AGAPstr] = lineList
##        tAGAP2 = time()
##        print '%s added to matrix: %s seconds.' % (AGAPstr,str(tAGAP2-tAGAP1))
    t2 = time()
    
    print 'readInMatrixFromFile() took %s seconds \n\tand created dict of length %s.\n' % (str(t2-t1), len(dictOfAGAPS))
    return [dictOfAGAPS, motifList]


def makeFwdAndRevCompRegExObj_IUPAC(motif, equals=0):
    import re
    
    motifPair = [motif, JamesDefs.revComp(motif)]
        
    targetContainsMotif  = '(%s|%s)' % (motifPair[0], motifPair[1])
    targetISMotif        = '^(%s|%s)&' % (motifPair[0], motifPair[1])
    
    if equals == 1:
        motif = targetISMotif
    elif equals == 0:
        motif = targetContainsMotif
    
    fwdRevComp_regExObj = re.compile(motif, re.IGNORECASE)
    return fwdRevComp_regExObj

def makeFwdAndRevCompRegExObj(motif):
    import re
    
    motifPair = [motif, JamesDefs.revComp(motif)]
    
    #  convert iupac string to regEx string
    for i in range(len(motifPair)):
        motifPair[i] = JamesDefs.iupac2regex(motifPair[i])
        
    motif = '(%s|%s)' % (motifPair[0], motifPair[1])
    
    fwdRevComp_regExObj = re.compile(motif, re.IGNORECASE)
    return fwdRevComp_regExObj



def findAllMotifsInAGAP(motifList, AGAP, dictOfFastas, matrixByAGAP):
        
    """
    populate the list of motif presensce/absence (1/0) values
    """
    listOfMotifPresenceVals = []
    
    for motif in motifList:
        #  compile regEx obj of motif to find both fwd or revComp versions of the motif
        motif = makeFwdAndRevCompRegExObj(motif)
        
        #  search AGAP and append 1 or 0 to list of motifs
        if motif.search(dictOfFastas[AGAP].seq.tostring()):
            listOfMotifPresenceVals.append(1)
        else:
            listOfMotifPresenceVals.append(0)
            
        matrixByAGAP[AGAP] = listOfMotifPresenceVals
    


    
   