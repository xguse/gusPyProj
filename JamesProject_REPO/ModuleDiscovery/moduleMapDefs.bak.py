#==============================================================
def countModuleInAll(module, dictWithModules):
    """
    Takes:	- Module (set of motifs)
                - dict of AGAPs with modules found as vals

    Does:	- resets each AGAPs moduleCount in dictWithModules[AGAPname][1] to 0
                - uses Set comparisons to loop over dictWithModules and update dictWithModules[AGAPname][1] in each AGAP that contains the module				
    
    Returns:	- totalHits and modifies dictWithModules[AGAPname][1] = 1 ,if it contains the current module, by reference
    """    
    

    totalHits = 0
    
    for AGAP in dictWithModules:
        #  initiate motifPresent attrib to 0
        dictWithModules[AGAP][1] = 0
        
        
        if module in dictWithModules[AGAP][0]:
            # update motifPresent attrib to reflect found motif
            dictWithModules[AGAP][1] = 1
            print 'Module = %s\nEntry for %s:\nModules: %s\nHits: %s\n' % (str(module), AGAP, str(dictWithModules[AGAP][0]), str(dictWithModules[AGAP][1]))
            totalHits += 1  
    
    return totalHits

    

#==============================================================
def buildModulesForAGAP(AGAPentryOfSortedMotifs, dictOfMotifSetsByAGAP, allMotifsEver, windowLen=500):
    '''
    takes    - AGAPdict entry
             - Dict to add list of motif sets(modules) to 
             - Length of widow (default to 500)
    
    does     - Creates a list of all motifs whose start pos is within windowLen, starting with smallest starting pos.
             - 
             - Convert that list into a Set obj which will allow comparison and also create nr version of the list for me and append that set to listOfMotifSets.
             - Repeat process with next smallest start pos.
            
    returns  - Nothing.  Simply adds entries in key[AGAPname]:val[listOfMotifSets] format to Dict included in call args.
    '''
    from sets import Set
    from probstat import Combination
    
    AGAPname = AGAPentryOfSortedMotifs[0]
    listOfSortedMotifs = AGAPentryOfSortedMotifs[1]
    listOfMotifSets = []
    leng = len(listOfSortedMotifs)
    motifIndex = 0
    while motifIndex < len(listOfSortedMotifs):
        motifsInWindow = []
        
        # add motif that starts window
        motifsInWindow.append(listOfSortedMotifs[motifIndex][0])
        
        # as long as the next motif's location minus the first motif's location is <= window length continue to append
        #     the next motif to motifsForModuleClass
        motifCounter = 1
        
        while motifIndex+motifCounter < len(listOfSortedMotifs)\
              and listOfSortedMotifs[motifIndex+motifCounter][1] - listOfSortedMotifs[motifIndex][1] < windowLen:
            listIndex = motifIndex+motifCounter
            motifsInWindow.append(listOfSortedMotifs[motifIndex+motifCounter][0])
            motifCounter += 1
        
        listOfMotifSets.append(Set(motifsInWindow))   
        motifIndex += 1
    
    # use Set() to get nr set of motif sets then iterate through the sets generating all posible combinations 
    # of all lengths in each set(since the motif sets have already passed the windowLen test, all combinations generated from these sets
    # must also be no more than windowLen bp apart) 
    setOfMotifSets = Set(listOfMotifSets)
    
    for motifSet in setOfMotifSets:
        listFromSet = list(motifSet)
        numOfMotifs = len(listFromSet)
        allCombinations = []
        comboLen = 1
        while comboLen <= numOfMotifs and comboLen <= 3:
            comboObj = Combination(listFromSet,comboLen)
            for item in comboObj:
                allCombinations.append(Set(item))
            comboLen+=1
        
  
        listOfMotifSets += allCombinations
        
        
    
    ##leng = len(listOfMotifSets)    
    listOfMotifSets = Set(listOfMotifSets)
    ##leng = len(listOfMotifSets) 
    #  assign listOfMotifSets to its AGAP key as part of a list with index[1] representing 
    #  the presense<1> or absence<1> of a module to be searched later
    dictOfMotifSetsByAGAP[AGAPname]= [list(listOfMotifSets), 0]
    #  add to list of all discovered motif sets then Set-ify and re-List-ify
    allMotifsEver+= list(listOfMotifSets)
    lenAll = len(allMotifsEver)
    w = 0
    
    for each in listOfMotifSets:
        print each
        
        
    print 'buildModulesForAGAP found %s modules for %s.' % (str(len(listOfMotifSets)), AGAPname)


#==============================================================

def spawnMotifInstances4AGAP(an_AGAPs_listOfMotifs):
    '''
    takes    - motifMap of one AGAP
    
    does     - creates instances(list of atribs [name,strand*,startPos]) of each occurence of each motif  *dont think i am gonna use strand
             - sorts the instances by startPos
            
    returns  - list: [AGAPname,[sortedListOfMotifInstances]]
    '''
    
    #  give list its AGAP name as index[0]
    AGAP_and_instances = [an_AGAPs_listOfMotifs[0][0]]
    instancesList = []
    while an_AGAPs_listOfMotifs:
        
        ## pop first two indexes of list (motif and motif_rc) to temp list
        motifsFwdRev = []
        motifsFwdRev.append(an_AGAPs_listOfMotifs.pop(0))
        motifsFwdRev.append(an_AGAPs_listOfMotifs.pop(0))
        
        ## next few steps: Combine fwd and rev lists of positions and check to see if the list is NON-empty
        motifName = motifsFwdRev[0][1]
                
        ##  Convert each pos into a instance list [motifName, pos] and add all to motifInstList only
        ##    if there are occurences of the motif
        if motifsFwdRev[0][2:]:
            for pos in motifsFwdRev[0][2:]:
                instancesList.append([motifName, int(pos)])
        if motifsFwdRev[1][2:]:
            for pos in motifsFwdRev[1][2:]:
                instancesList.append([motifName, int(pos)])
        
    #  Sort instanceList append it to AGAP_and_instances and return that list
    instancesList.sort(cmp=lambda x,y: x[1]-y[1])    
    AGAP_and_instances.append(instancesList)
    
    return AGAP_and_instances
    
    