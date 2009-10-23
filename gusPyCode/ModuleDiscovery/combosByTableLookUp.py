from gusPyCode.defs.defs_moduleByTableLookUp import *
from gusPyCode.defs.JamesDefs import groupByField
from gusPyCode.defs.doug_hypergeometric import hyperGeoPvalue
from time import time
from sets import Set


#========================= User Defined Variables =========================
matrixFile      = '/Users/biggus/Documents/MBGB/Rotations/James/Data/MotifMaps/matrixOfMotif_vs_AGAP.txt'
clusterDefs     = '/Users/biggus/Documents/MBGB/Rotations/James/Data/ClusterDefs/053008_incompleteTimeCourse.txt'
motifList       = '/Users/biggus/Documents/MBGB/Rotations/James/Data/DegenMotifs/aedesAnopheles7mer.TSS.nr.motif' #'/Users/biggus/Documents/MBGB/Rotations/James/Data/testData4Mapper/5motifsForMapping.motif'  

outFile = open('/Users/biggus/Documents/MBGB/Rotations/James/Data/MotifMaps/motifCombos.txt','w')

#==========================================================================

#  Reconsitute the matrix into memory
matrixAndList = readInMatrixFromFile(matrixFile)
matrixOfAGAPsvMotifs = matrixAndList[0]

#  Initialize and group clusterDef list 
clusterDefs = map(lambda line : line.strip(), open(clusterDefs, 'rU').readlines())
clusterDefs = groupByField(clusterDefs,0)

#  Initialize motifList
motifList = map(lambda line : line.strip(), open(motifList, 'rU').readlines())

#  Motif combos
motifCombos = generateMotifCombos(motifList,2)

#  Counting dictionary to hold presence/absence for each combo in AGAP
countingDict = {}

for AGAP in matrixOfAGAPsvMotifs:
    countingDict[AGAP] = 0

#  Motif-2-matrixIndex translation dict
motif2indexDict = {}
i=0
for motif in motifList:
    motif2indexDict[motif] = i
    i+=1

#  To catch params for AGAP/combo    
paramsForHyperGeo = []

#  big N = Pop
N = len(matrixOfAGAPsvMotifs)

#  to catch any AGAPs that happen to be in the cluster Defs but not in our list of boundry regions
AGAPsNotFound = Set()

# Loop through 
for combo in motifCombos:
    t1 = time()
    #  Clear all AGAP's slots in counting dict to ensure no leftover counts
    for AGAP in countingDict:
        countingDict[AGAP] = 0
    
    #  Count combo in AGAP
    for AGAP in matrixOfAGAPsvMotifs:
        findMotifComboInAGAP(AGAP, matrixOfAGAPsvMotifs, combo, countingDict, motif2indexDict)
        
    #  big K = how many hits in Pop
    K = getBigK(countingDict)
    
    for cluster in clusterDefs:
        
        #  little n = size of cluster
        n = len(cluster)
        
        #  little k = how many hits in cluster
        k = getLittleK(countingDict, cluster, AGAPsNotFound)
        
        #  format params into line for paramsForHyperGeo -----> 'combo;clusterName    N    n    K    k'
        params = '%s;%s\t%s\t%s\t%s\t%s' % ('_'.join(combo),cluster[0][0],str(N),str(n),str(K),str(k))
        paramsForHyperGeo.append(params)
    
    t2 = time()
    print "counting for %s took %s sec." % (str(combo), str(t2-t1))
    
print 'Durring the counting I encountered the following AGAPs that were\nin the Cluster Defs but not in the matrix:\n%s' \
      % '\n'.join(list(AGAPsNotFound))

outFile.write('\n'.join(paramsForHyperGeo))