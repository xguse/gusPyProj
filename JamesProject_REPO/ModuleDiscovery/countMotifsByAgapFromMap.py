from sets import Set
from time import time
from JamesDefs import *
from moduleMapDefs import *
from probstat import Combination

"""
takes:    - file with map of given motif locations
          - cluster definition file
          - path to out
does:
returns:
"""


#--------- Script Specific Function Definitions ---------------------

#--------------------------------------------------------------------



#========================= User Defined Variables =========================

clusterDefPath             = '/Users/biggus/Documents/MBGB/Rotations/James/Data/ClusterDefs/052708_3xTukeyHSD_tissueClusters.txt'
motifMapPath               = '/Users/biggus/Documents/MBGB/Rotations/James/Data/testData4Mapper/MGandUpAt24.map'
nrDegenMotifsFromMDOS_path = '/Users/biggus/Documents/MBGB/Rotations/James/Data/DegenMotifs/aedesAnopheles7mer.TSS.nr.motif'
outFile                    = '/Users/biggus/Documents/MBGB/Rotations/James/Data/testData4Mapper/MGandUpAt24.matrix'



#==========================================================================

#  Initialize and format clusterDef list of lists by clusterName
clusterDefList = map(lambda line: line.strip(), open(clusterDefPath, 'rU').readlines())
clusterDefList = groupByField(clusterDefList, 0)

#  Initialize and format list of lists by AGAP
motifMap = map(lambda line: line.strip(), open(motifMapPath, 'rU').readlines())
print 'motifMap has %s items.' % (len(motifMap))

startMotifMapTime = time()
motifMap = groupByField_silent(motifMap, 0)      
endMotifMapTime = time()
print 'groupByField on motifMap took %s min.' % (str((endMotifMapTime-startMotifMapTime)/60))

#  Initialize MDOS_degenMotifs
MDOS_degenMotifs = map(lambda line: line.strip(), open(nrDegenMotifsFromMDOS_path, 'rU').readlines())

#  Convert motif map into Dict[AGAPs] of lists['instances of motifs(each should be list of IUPAC_motifStr;startPos)']
#  Then collapse this to be simply a nonRedundant list of each motif found in this AGAP
map2DictStrt = time()
dictOfMotifsByAGAP = {}

for AGAPList in motifMap:
    oneAGAPsMotifs = spawnMotifInstances4AGAP(AGAPList)
    nrListOfMotifs = Set()
    for motifInstance in oneAGAPsMotifs[1]:
        nrListOfMotifs.add(motifInstance[0])
    oneAGAPsMotifs[1] = nrListOfMotifs
    dictOfMotifsByAGAP[oneAGAPsMotifs[0]] = oneAGAPsMotifs[1]
    ##len_dictOfMotifsByAGAP = len(dictOfMotifsByAGAP)
map2DictEnd = time()    
print 'map2dict took %s min.' % (str((map2DictEnd-map2DictStrt)/60))


#  Create dict to catch list of motifs and whether they exist in the AGAP
dictOfMotifPresenceByAGAP = {}

#  Loop to check for motif in each AGAP and populate list = [[Motif1, 0|1], [Motif2, 0|1]] / 0 is no, 1 is yes
for AGAP in dictOfMotifsByAGAP:
   
    #  Define each AGAP's val to be a list so that we can append the motif to it later
    dictOfMotifPresenceByAGAP[AGAP] = []
    for motif in MDOS_degenMotifs:

        if motif in dictOfMotifsByAGAP[AGAP]:
            
            dictOfMotifPresenceByAGAP[AGAP].append([motif, 1])
        
        else:
            
            dictOfMotifPresenceByAGAP[AGAP].append([motif, 0])



#  Construct a list representing a tab delimited output which will be written to a file later

outPutList = []

firstLine = '\t'.join(MDOS_degenMotifs)
firstLine = '\t'+firstLine
outPutList.append(firstLine+'\n')

for AGAP in dictOfMotifPresenceByAGAP:
    AGAP_line = AGAP
    for motif in dictOfMotifPresenceByAGAP[AGAP]:
        AGAP_line+= '\t'+str(motif[1])
    
    outPutList.append(AGAP_line+'\n')


outFile = open(outFile, 'w').writelines(outPutList)
   
print 'yay'

