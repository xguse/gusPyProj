from sets import Set
from time import time
from gusPyCode.defs.JamesDefs import *
from gusPyCode.ModuleDiscovery.moduleMapDefs import *
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
outFile                    = 'path/To/out'

mostMotifsInSet = 3

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
map2DictStrt = time()
dictOfSortedMotifsByAGAP = {}
for AGAPList in motifMap:
    oneAGAPsMotifs = spawnMotifInstances4AGAP(AGAPList)
    dictOfSortedMotifsByAGAP[oneAGAPsMotifs[0]] = oneAGAPsMotifs[1]
    len_dictOfSortedMotifsByAGAP = len(dictOfSortedMotifsByAGAP)
map2DictEnd = time()    
print 'map2dict took %s min.' % (str((map2DictEnd-map2DictStrt)/60))

#  Generate list of all motif combinations(modules) of lengths between 1 and mostMotifsInSet
listOfSearchModules = []
comboLen = 1
lenOfDegenMotifs = len(MDOS_degenMotifs)
tl1 = time()
while comboLen <= lenOfDegenMotifs and comboLen <= mostMotifsInSet:
    comboObj = Combination(MDOS_degenMotifs,comboLen)
    for item in comboObj:
        listOfSearchModules.append(Set(item))
    comboLen+=1

len_listOfSearchModules = len(listOfSearchModules)
tl2 = time()
x = 0
print "Completed listOfSearchModules: %s min" % (str((tl2-tl1)/60))

#  Create dict to catch list of modules(modules = Set([listOfMotifsInWindowSize])) 
#  Set up loop through dictOfSortedMotifsByAGAP to send each AGAP rec in turn to buildModulesForAGAP()
#  buildModulesForAGAP() will loop through list of motifs to build modules and add the list of sets to new dict using AGAPname as key
#  ***default windowLen = 500 for buildModulesForAGAP() -> its an optional arg and changeable***
dictOfMotifSetsByAGAP = {}

t1 = time()
for k,v in dictOfSortedMotifsByAGAP.iteritems():
    buildModulesForAGAP([k,v], dictOfMotifSetsByAGAP)

t2 = time()

print "dictOfSortedMotifsByAGAP loop took %s minutes." % (str((t2-t1)/60))

#  Count 

    
print 'yay'

