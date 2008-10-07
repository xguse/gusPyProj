from sets import Set
from time import time
from JamesDefs import *
from moduleMapDefs import *
from probstat import Combination
import sys
# py modules under development

"""
takes:    - file with map of given motif locations
          - cluster definition file
          - path to out
does:
returns:
"""

#========================= User Defined Variables =========================

clusterDefPath             = '/Users/biggus/Documents/James/Data/ClusterDefs/7-08_Descriptive_time-course_Clusters.txt'
motifMapPath               = '/Users/biggus/Documents/James/Data/2KB/2kb_Sequence/2kb_Combo/2Kb_AllMosquitoes/MosqMotifs/upstream_exclsv-conserved_mosquito-motifs_nr.smLine.map'
nrDegenMotifsFromMDOS_path = '/Users/biggus/Documents/James/Data/2KB/2kb_Sequence/2kb_Combo/2Kb_AllMosquitoes/MosqMotifs/upstream_exclsv-conserved_mosquito-motifs_nr.txt'
dictOfMotifSetsByAGAP_OUT  = open('/Users/biggus/Documents/James/Data/2KB/2kb_Sequence/2kb_Combo/2Kb_AllMosquitoes/fake1.txt','w')
outFile                    = open('/Users/biggus/Documents/James/Data/2KB/2kb_Sequence/2kb_Combo/2Kb_AllMosquitoes/fake2','w')


mostMotifsInSet = 2

#==========================================================================
#--------- Script Specific Function Definitions ---------------------
def buildMotifPairsForAGAP(AGAPentryOfSortedMotifs, windowLen=500):
    '''
    takes    - AGAPdict entry
             - Length of widow (default to 500)
    
    does     - calculates nr list of all motif pairs that are less than <windowLen> apart
             - Current version uses 'Sets' to remove redundancy
  
    returns  - List in form: [AGAPname, [listOfMotifSets,0]].  The '0' acts as a presence/absence switch later on
    '''    
    from sets import Set,ImmutableSet

    AGAPname = AGAPentryOfSortedMotifs[0]
    listOfSortedMotifs = AGAPentryOfSortedMotifs[1]
    len_listOfSortedMotifs = len(listOfSortedMotifs)
    setOfMotifSets = Set()

    #  	Gonna try a 'conveyor-belt' type method of popping the first motif and compiling a list of those that satisfy the tests.
    #	That way I should not have to use complicated methods of keeping track of where I am in the stream    
    motifPairList = []
    
    while listOfSortedMotifs:
        primaryMotif = listOfSortedMotifs.pop(0)
        
        #  hold a copy of all motif entries satisfying the window size requirement
        motifsInWindow = []
        
        # copy motifs in window to working list
        # BREAK once window is exited to avoid extra computations
        for each in listOfSortedMotifs:
            if each[1] <= primaryMotif[1]+windowLen:
                motifsInWindow.append(each)
            else:
                break
            
        # after window-test, test for overlap
        # I will simply (1) append each succesful test to the end of the list
        #    and then (2) delete th original indexes, leaving only those that passed
        len_mIW = len(motifsInWindow)
        for i in range(0,len_mIW):
            
            nMotifStart        = motifsInWindow[i][1]
            len_primaryMotif   = len(primaryMotif[0])
            primaryMotifStart  = primaryMotif[1]
            result             = nMotifStart-len_primaryMotif
            
            if result > primaryMotifStart:
                motifsInWindow.append(motifsInWindow[i])       
        del motifsInWindow[0:len_mIW]  #  (Step2) from comments above
        
        # Add motif sets with primaryMotif to SetofMotifSets
        for m in motifsInWindow:
            setOfMotifSets.add(ImmutableSet([primaryMotif[0],m[0]]))
            
    if len(setOfMotifSets) > 0:
        #print 'buildMotifPairsForAGAP found %i pairs for %s.' % (len(setOfMotifSets), AGAPname)
        return [AGAPname,[list(setOfMotifSets),0]]
    
    #  cant list(None) or len(None) so I must manually return an empty list and set len to zero
    else:   
        print 'buildMotifPairsForAGAP found %i pairs for %s.' % (0, AGAPname)
        return [AGAPname,[[],0]]
        


#-=-=-=-=-=-=-=-=-=-=-=-=-=-


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
    
def groupMap(listOfTabbedStrings, fieldGroupedBy): 
    """ WARNING!! listOfTabbedStrings will be destroyed!!
    Example:
    fieldGroupedBy = 0
    listOfTabbedStrings = ['a\t1','a\t2','b\t1']
    result-> [[['a','1'],['a','2']],[['b','1']] ]
    """

    returnList = []
    #  allows us to test current group class
    currentGroup = listOfTabbedStrings[0].split('\t')
    
    tempList = []
    i = 0
    while i < len(listOfTabbedStrings):
        
        
        if listOfTabbedStrings[i].split('\t')[0] == currentGroup[0]:
            tempList.append(listOfTabbedStrings[i].rstrip('\n').split('\t'))
            i+=1
        else:
            #print len(tempList)
            returnList.append(tempList)
            tempList = []
            currentGroup = listOfTabbedStrings[i].split('\t')
    return returnList
#--------------------------------------------------------------------

totalTimeStart = time()



##  Initialize and format clusterDef list of lists by clusterName
clusterDefList = map(lambda line: line.strip(), open(clusterDefPath, 'rU').readlines())
clusterDefList = groupByField(clusterDefList, 0)

#  Initialize and format list of lists by AGAP
motifMap = map(lambda line: line.strip(), open(motifMapPath, 'rU').readlines())
print 'motifMap has %s items.' % (len(motifMap))

startMotifMapTime = time()
cpu1 = cpu()
motifMap = groupMap(motifMap, 0)      
cpu2 = cpu()
endMotifMapTime = time()
print 'groupMap on motifMap produced %i groups and took %.3f min.' % (len(motifMap),(endMotifMapTime-startMotifMapTime)/60.0)

x = motifMap[0]

## To kill after time report for debugging
#sys.exit()

#  Initialize MDOS_degenMotifs
MDOS_degenMotifs = map(lambda line: line.strip(), open(nrDegenMotifsFromMDOS_path, 'rU').readlines())

#  Convert motif map into Dict[AGAPs] of lists['instances of motifs(each should be list of IUPAC_motifStr;startPos)']
map2DictStrt = time()
dictOfSortedMotifsByAGAP = {}
for AGAPList in motifMap:
    oneAGAPsMotifs = spawnMotifInstances4AGAPv2(AGAPList)
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
x = listOfSearchModules[0]


    

    
print "Completed listOfSearchModules: %.3f min.  (%i were generated)" % ((tl2-tl1)/60.0,len_listOfSearchModules)

#  Create dict to catch list of modules(modules = Set([listOfMotifsInWindowSize])) 
#  Set up loop through dictOfSortedMotifsByAGAP to send each AGAP rec in turn to buildModulesForAGAP()
#  buildModulesForAGAP() will loop through list of motifs to build modules and add the list of sets to new dict using AGAPname as key
#  ***default windowLen = 500 for buildModulesForAGAP() -> its an optional arg and changeable***
dictOfMotifSetsByAGAP = {}

t1 = time()
for k,v in dictOfSortedMotifsByAGAP.iteritems():
    ##dev_buildMotifPairsForAGAP([k,v], dictOfMotifSetsByAGAP)
    AGAPsResultList = buildMotifPairsForAGAP([k,v])
    if AGAPsResultList:
        if dictOfMotifSetsByAGAP.has_key(AGAPsResultList[0]):
            print "Warning! key:%s already present in 'dictOfMotifSetsByAGAP'!  Quiting."
            sys.exit()
        dictOfMotifSetsByAGAP[AGAPsResultList[0]] = AGAPsResultList[1]
    else:
        print "Warning! buildMotifPairsForAGAP for %s did not return a value! Quiting." % (k)
        sys.exit()

t2 = time()

print "dictOfSortedMotifsByAGAP loop took %s minutes." % ((t2-t1)/60.0)


writeOutDict(dictOfMotifSetsByAGAP,dictOfMotifSetsByAGAP_OUT)
x=1

# Delete dictOfSortedMotifsByAGAP to free memory:
del dictOfSortedMotifsByAGAP


#  list to hold ModuleClass;ClusterID;hyperGeoParams
hyperGeoParams_4_moduleClusterPairs = []

#  list to hold error output from 'for cluster in clusterDefList:' loop
notIn_dictOfMotifSetsByAGAP = []

m=0
for module in listOfSearchModules:
    modTime1 = time()
    modCpuTime1 = cpu()
    m+=1
    
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
        outFile.write(str(module)+','+seqID[0]+tab+str(len(dictOfMotifSetsByAGAP))+tab+str(len(cluster))+tab+str(moduleCountInAll)+tab+str(moduleCountInCluster)+'\n') 
        #outFile.flush()
    
    modTime2 = time()
    modCpuTime2 = cpu()
    print 'Module %i: clock = %.6f cpu = %.6f' % (m, modTime2-modTime1, modCpuTime2-modCpuTime1)
        
t2 = time()        

notIn_dictOfMotifSetsByAGAP = Set(notIn_dictOfMotifSetsByAGAP)
notIn_dictOfMotifSetsByAGAP = list(notIn_dictOfMotifSetsByAGAP)
#  report about missing AGAPs from dictOfMotifSetsByAGAP
print 'The following %s AGAPs were not found in dictOfMotifSetsByAGAP:' % (str(len(notIn_dictOfMotifSetsByAGAP)))
for each in notIn_dictOfMotifSetsByAGAP:
    print each






## Now writing as I go instead of:  outFile.writelines(hyperGeoParams_4_moduleClusterPairs)
outFile.close()

totalTimeEnd = time()

print 'The cluster counting took %.3f minutes.\nTotalTime: %.3f' \
      % ((t2-t1)/60.0, (totalTimeEnd-totalTimeStart)/60.0)

sys.exit('proper end')    


