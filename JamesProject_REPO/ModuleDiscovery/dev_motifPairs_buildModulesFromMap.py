from sets import Set
from time import time
from JamesDefs import *
from moduleMapDefs import *
from probstat import Combination
import sys
# py modules under development
from dev_buildMotifPairsForAGAP import buildMotifPairsForAGAP, dev_buildMotifPairsForAGAP

"""
takes:    - file with map of given motif locations
          - cluster definition file
          - path to out
does:
returns:
"""


#--------- Script Specific Function Definitions ---------------------
def writeOutDict(dict,fileHandle):
    outList = []
    print "Dict Size: %i" % (len(dict))
    for k1,v1 in dict.iteritems():
        recStr = k1
        for each in v1[0]:
            l = list(each)
            if len(l)>1:
                setStr = '\t%s:%s' % (l[0],l[1])
                recStr+=setStr
        recStr+='\n'
        outList.append(recStr)

    print 'ModulesByAGAP written: %i' % (len(outList))
    fileHandle.writelines(outList)
    x=1
#--------------------------------------------------------------------

totalTimeStart = time()

#========================= User Defined Variables =========================

clusterDefPath             = '/Users/biggus/Documents/James/Data/ClusterDefs/VERIFIED_FINAL_tissue-enriched.txt'
motifMapPath               = '/Users/biggus/Documents/James/Data/2KB/2kb_Sequence/2kb_Combo/2Kb_AllMosquitoes/MosqMotifs/collapsed_upstream_ex-conserved_mosquito-motifs_nr.copRvCmp.map'
nrDegenMotifsFromMDOS_path = '/Users/biggus/Documents/James/Data/2KB/2kb_Sequence/2kb_Combo/2Kb_AllMosquitoes/MosqMotifs/collapsed_upstream_ex-conserved_mosquito-motifs_nr.copRvCmp.txt'
dictOfMotifSetsByAGAP_OUT  = open('/Users/biggus/Documents/James/Data/2KB/2kb_Sequence/2kb_Combo/2Kb_AllMosquitoes/dictOfMotifPairsByAGAP_tissueEnriched.txt','w')
outFile                    = open('/Users/biggus/Documents/James/Data/2KB/2kb_Sequence/2kb_Combo/2Kb_AllMosquitoes/2KBup_collapsed_upstream_ex-conserved_mosquito-motifs_nr.te.hgp.txt','w')

mostMotifsInSet = 2

#==========================================================================

##  Initialize and format clusterDef list of lists by clusterName
clusterDefList = map(lambda line: line.strip(), open(clusterDefPath, 'rU').readlines())
clusterDefList = groupByField(clusterDefList, 0)

#  Initialize and format list of lists by AGAP
motifMap = map(lambda line: line.strip(), open(motifMapPath, 'rU').readlines())
print 'motifMap has %s items.' % (len(motifMap))
startMotifMapTime = time()
motifMap = groupByField_silent(motifMap, 0)      
endMotifMapTime = time()
print 'groupByField on motifMap took %.3f min.' % ((endMotifMapTime-startMotifMapTime)/60.0)

## To kill after time report for debugging
#sys.exit()

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
print 'map2dict took %.3f min.' % ((map2DictEnd-map2DictStrt)/60.0)

#  Generate list of all motif combinations(modules) of lengths between 1 and the user defined mostMotifsInSet
listOfSearchModules = []
comboLen = 2
lenOfDegenMotifs = len(MDOS_degenMotifs)
tl1 = time()

comboObj = Combination(MDOS_degenMotifs,comboLen)
for item in comboObj:
    listOfSearchModules.append(Set(item))


len_listOfSearchModules = len(listOfSearchModules)
tl2 = time()
x = 0
print "Completed listOfSearchModules: %.3f min.  (%i were generated)" % ((tl2-tl1)/60.0,len_listOfSearchModules)

#  Create dict to catch list of modules(modules = Set([listOfMotifsInWindowSize])) 
#  Set up loop through dictOfSortedMotifsByAGAP to send each AGAP rec in turn to buildModulesForAGAP()
#  buildModulesForAGAP() will loop through list of motifs to build modules and add the list of sets to new dict using AGAPname as key
#  ***default windowLen = 500 for buildModulesForAGAP() -> its an optional arg and changeable***
dictOfMotifSetsByAGAP = {}

t1 = time()
for k,v in dictOfSortedMotifsByAGAP.iteritems():
    dev_buildMotifPairsForAGAP([k,v], dictOfMotifSetsByAGAP)
    ##AGAPsResultList = buildMotifPairsForAGAP([k,v])
    ##if AGAPsResultList != None:
        ##if dictOfMotifSetsByAGAP.has_key(AGAPsResultList[0]):
            ##print "Warning! key:%s already present in 'dictOfMotifSetsByAGAP'!  Quiting."
            ##sys.exit()
        ##dictOfMotifSetsByAGAP[AGAPsResultList[0]] = AGAPsResultList[1]

t2 = time()

print "dictOfSortedMotifsByAGAP loop took %s minutes." % ((t2-t1)/60.0)


writeOutDict(dictOfMotifSetsByAGAP,dictOfMotifSetsByAGAP_OUT)
x=1


#  list to hold ModuleClass;ClusterID;hyperGeoParams
hyperGeoParams_4_moduleClusterPairs = []

#  list to hold error output from 'for cluster in clusterDefList:' loop
notIn_dictOfMotifSetsByAGAP = []

m=0
for module in listOfSearchModules:
    m+=1
    print 'Module '+str(m)
    #  Count how many AGAPs in total list have module 
    moduleCountInAll = None
    moduleCountInAll = countModuleInAll(module, dictOfMotifSetsByAGAP)
 
    
    for cluster in clusterDefList:

        #  count how many seq in Cluster have motifPresent set to 1
        moduleCountInCluster = 0
        
        for seqID in cluster:
            if dictOfMotifSetsByAGAP.has_key(seqID[1]):
                moduleCountInCluster += dictOfMotifSetsByAGAP[seqID[1]][1]
            else:
                notIn_dictOfMotifSetsByAGAP.append(seqID[1]) #+' not present in dictOfMotifSetsByAGAP dictionary!')
            
        #  Format hyperGeoParams_4_moduleClusterPairs record
        tab = '\t'
        # ------------------------------------->  N,n,K,k
        hyperGeoParams_4_moduleClusterPairs.append(str(module)+','+seqID[0]+tab+str(len(dictOfMotifSetsByAGAP))+tab+str(len(cluster))+tab+str(moduleCountInAll)+tab+str(moduleCountInCluster)+'\n') 

t2 = time()        

notIn_dictOfMotifSetsByAGAP = Set(notIn_dictOfMotifSetsByAGAP)
notIn_dictOfMotifSetsByAGAP = list(notIn_dictOfMotifSetsByAGAP)
#  report about missing AGAPs from dictOfMotifSetsByAGAP
print 'The following %s AGAPs were not found in dictOfMotifSetsByAGAP:' % (str(len(notIn_dictOfMotifSetsByAGAP)))
for each in notIn_dictOfMotifSetsByAGAP:
    print each


t3 = time()
#c = 0
#while c < len(hyperGeoParams_4_moduleClusterPairs):
    ##  explode string to list of params
    #params = hyperGeoParams_4_moduleClusterPairs[c].split('\t')
    
    ##  calc hyperGeo pVal
    #pVal = hyperGeoPvalue(int(params[1]),int(params[2]),int(params[3]),int(params[4]))
    #print pVal
    #hyperGeoParams_4_moduleClusterPairs[c] = hyperGeoParams_4_motifClusterPairs[c]+'\t'+str(pVal)+'\n'
    #c+=1    
t4 = time()



outFile.writelines(hyperGeoParams_4_moduleClusterPairs)
outFile.close()

totalTimeEnd = time()

print 'The cluster counting took %.3f minutes.\nThe hyperGeo function took %.3f minutes.\nTotalTime: %.3f' \
      % ((t2-t1)/60.0, (t4-t3)/60.0, (totalTimeEnd-totalTimeStart)/60.0)

sys.exit('proper end')    


