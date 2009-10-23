import re
from MDOSX_classes import OrthoGroup
from TAMO.seq import Fasta
from TAMO import MotifMetrics

def searchOrthoGroups(listOfMotifObjs, dictOfOrthoGroups, scoreFactor=None):
    """Takes listOfMotifObjs<list> and dictOfOrthoGroups<dict> and feeds each motifObj
    to to each orthoGroupObj's searchMotif method."""
    
    # validation
    ## to be added
    
    if scoreFactor == None:
        scoreFactor = 0.75
    
    for motif in listOfMotifObjs:
        for k in dictOfOrthoGroups:
            dictOfOrthoGroups[k].searchMotif(motif,scoreFactor)



# ######################
# ######################


def spawnOrthoGroups(promoterFileList,nWayOrthoList):
    """Takes promoterFileList<listOfPaths> and nWayOrthoList<listOfLists> and spawns the orthoGroup
    objects in a dictionary with keys = 'geneName1:geneName2:etc' that will be used to run the combined
    hypergeometric analysis."""
    
    
    
    # validation
    assert type(promoterFileList) == type([]), \
           '''promoterFileList must be a list of file paths.
           You provided type: "%s"'''\
           % (type(promoterFileList))
    assert type(promoterFileList[0]) == type(''), \
           '''promoterFileList must be a list of file paths.
           promoterFileList[0] != type(''): "%s"'''\
           % (type(promoterFileList[0]))
    
    # load promoters
    allPromoters = {}
    for i in range(len(promoterFileList)):
        oneGenome = Fasta.file2dict(promoterFileList[i])
        for j in oneGenome:
            allKeys = allPromoters.keys()
            assert j not in allKeys, \
                   '''Detected duplicate gene name in promoterFileList! "%s"'''\
                   % (j)
            allPromoters[j] = oneGenome[j]
    
    # Build Groups
    orthoGroups = {}
    for i in range(len(nWayOrthoList)):
        groupDict = {}
        for j in range(len(nWayOrthoList[i])):
            if allPromoters[nWayOrthoList[i][j]]:
                groupDict[nWayOrthoList[i][j]] = allPromoters[nWayOrthoList[i][j]]
            else:
                break # we do not want orthoGroups that are missing members
        
        if len(groupDict) != len(nWayOrthoList[i]):
            break # we do not want orthoGroups that are missing members
        else:
            nWayOrthoList[i].sort()
            orthoGroups[':'.join(nWayOrthoList[i])] = OrthoGroup(groupDict)
            
    return orthoGroups


# ######################
# ######################

def intxnOf_N_Sets(listOfSets):
    intxn = None
    for i in range(len(listOfSets)):
        if i+1 == len(listOfSets):
            return intxn
        else:
            intxn = listOfSets[i]&listOfSets[i+1]

# ######################
# ######################        

def parseCoRegdGeneList(coRegdGeneListFile):
    # read coRegdGeneListFile
    coRegdGeneListFile = map(lambda l: l.strip(), open(coRegdGeneListFile,'r').readlines())
    # remove blank lines
    coRegdGeneList = []
    for l in coRegdGeneListFile:
        if l != '':
            coRegdGeneList.append(l)
    
    assert coRegdGeneList[0].startswith('>'), \
           """coRegdGeneListFile is expected to be a fasta-like file starting with '>GenomeName'.
           Your first line looks like this: '%s'""" % (coRegdGeneList[0])
    
    genesByGenome = map(lambda l: l.rstrip().split('\n'), '\n'.join(coRegdGeneList).split('>'))
    genesByGenome.pop(0)
    for i in genesByGenome:
        i.pop(0)
    
    return genesByGenome
            
# ######################
# ######################

def buildCoRegdOrthoGroups(nWayOrthoList, coRegdGeneListFile):
    """\tTakes nWayOrthoList<list> and coRegdGeneListFile<'filePath'> and returns the subset of
    orthoGroups in nWayOrthoList that share the expression profile of interest."""
    
    genesByGenome = parseCoRegdGeneList(coRegdGeneListFile)
    nWayOrthoList = nWayOrthoList[:] # ensure I dont alter the list outside of the def
    
    # Replace gene names in each genome list that have n-way orthologs with orthoGroup
    for i in range(len(genesByGenome)):
        for j in range(len(genesByGenome[i])):
            for nWayOrtho in nWayOrthoList:
                if genesByGenome[i][j] in nWayOrtho:
                    genesByGenome[i][j] = nWayOrtho[:]
                    continue
    # Then taking the intersection of all lists gives you your coRegdOrthos
    # set-ify all geneNames or OrthoGroups in all genomeLists
    for i in range(len(genesByGenome)):
        for j in range(len(genesByGenome[i])):
            genesByGenome[i][j] = frozenset(genesByGenome[i][j])
    # set-ify the genomeLists
    for i in range(len(genesByGenome)):
        genesByGenome[i] = frozenset(genesByGenome[i])
    
    intxn = intxnOf_N_Sets(genesByGenome)
    
    # Convert back into a list of alpha-sorted lists
    intxn = list(intxn)
    for i in range(len(intxn)):
        intxn[i] = list(intxn[i])
        intxn[i].sort()
        
    return intxn
    

    


def combineOrthoTabs(listOfOrthoTabFiles,numOfSpecies):
    """\tTake list of files containing ortholog tables in format:GeneIDspecies1[tab]GeneIDspecies2
    Return a list of orthologs found in all species."""
    
    orthoTabs     = []
    for i in listOfOrthoTabFiles:
        data = map(lambda l: l.strip().split('\t'), open(i,'rU').readlines())
        orthoTabs.append(data[:])
    
    # Move longest list to first position and use _IT_ as root list for merging
    orthoTabSizes = []
    for i in orthoTabs:
        orthoTabSizes.append(len(i))
    ##print '%s' % (str(orthoTabSizes))
    for i in range(len(orthoTabs)):
        if len(orthoTabs[i]) == max(orthoTabSizes):
            biggest = orthoTabs.pop(i)
            orthoTabs.insert(0,biggest)
    
    orthoTabSizes2 = []
    for i in orthoTabs:
        orthoTabSizes2.append(len(i))       
    ##print '%s' % (str(orthoTabSizes2))
    
    # Set-ify all orthoPairs
    for i in range(len(orthoTabs)):
        for j in range(len(orthoTabs[i])):
            orthoTabs[i][j] = set(orthoTabs[i][j])
    
    # Merge sets if any geneIDs are in common
    for i in range(1,len(orthoTabs)):                          # Start w/ 2nd orthoList and compare it and all after to 1st list
        for j in range(len(orthoTabs[i])):                     # For every j in 2nd or higher 
            for k in range(len(orthoTabs[0])):                 # Test every k in 1st orthoList for j
                if orthoTabs[0][k]&orthoTabs[i][j]:            # If there is intersection bt k and ij,
                    orthoTabs[0][k].update(orthoTabs[i][j])    # merge k and ij
                    
    # Remove redundancy and rename merged list for ease of use
    print 'len of mergedIDs list = %s' % (len(orthoTabs[0]))
    mergedIDs = orthoTabs[0]
    for i in range(len(mergedIDs)):
        mergedIDs[i] = frozenset(mergedIDs[i])
    mergedIDs = list(set(mergedIDs))
    print 'len of nr mergedIDs list = %s' % (len(mergedIDs))
    
    
    # Sort each orthoSet
    for i in range(len(mergedIDs)):
        mergedIDs[i] = list(mergedIDs[i])
        mergedIDs[i].sort()
        
    # Retain only those relationships that are the correct length
    # AND do not repeat Species prefixs ('AGAP')
    filteredIDs = []
    for orthoSet in mergedIDs:
        prefixes = []
        for ID in orthoSet:
            prefixes.append(ID[:4])
        if len(prefixes) != numOfSpecies:
            continue
        if len(prefixes) == len(set(prefixes)):
            filteredIDs.append(orthoSet)
            
 
            
    return filteredIDs
            
changeLog = """2009-03-29 -- started log (combineOrthoTabs was already written)
2009-03-29 -- added buildCoRegdOrthoGroups() [complete]
2009-03-29 -- added parseCoRegdGenes() [complete]
2009-03-29 -- added intxnOf_N_Sets() [complete]
2009-03-30 -- added spawnOrthoSets() [functional -> add info-logging?]
2009-03-30 -- added searchOrthoGroups() [needs validation code]
"""
