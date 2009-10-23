#!/usr/bin/env python



def buildMotifPairsForAGAP(AGAPentryOfSortedMotifs, windowLen=500):
    '''
    takes    - AGAPdict entry
             - Length of widow (default to 500)
    
    does     - calculates nr list of all motif pairs that are less than <windowLen> apart
             - Current version uses 'Sets' to remove redundancy
  
    returns  - List -> [AGAPname, [listOfMotifSets]].
    '''    
    from sets import Set

    AGAPname = AGAPentryOfSortedMotifs[0]
    listOfSortedMotifs = AGAPentryOfSortedMotifs[1]
    len_listOfSortedMotifs = len(listOfSortedMotifs)
    ListOfMotifPairSets = []

    #  	Gonna try a 'conveyor-belt' type method of popping the first motif and compiling a list of those that satisfy the tests.
    #	That way I should not have to use complicated methods of keeping track of where I am in the stream
    
    motifPairList = []
    
    while listOfSortedMotifs:
        primaryMotif = listOfSortedMotifs.pop(0)
        
        #  hold a copy of all motif entries satisfying the window size requirement
        motifsInWindow = []
        
        # copy motifs in window to working list
        for each in listOfSortedMotifs:
            if each[1] <= primaryMotif[1]+500:
                motifsInWindow.append(each)



























def buildModulesForAGAPv2(AGAPentryOfSortedMotifs, dictOfMotifSetsByAGAP, windowLen=500):
    '''
    takes    - AGAPdict entry
             - Dict to add list of motif sets(modules) to 
             - Length of widow (default to 500)
    
    does     - NOTE: need to revise this b/c i have made changes (current max combo len = 4 for example)
             - Creates a list of all motifs whose start pos is within windowLen, starting with smallest starting pos.
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
   
    #  this is just to let me look at the list in term
    #for mtf in listOfSortedMotifs:
        #print mtf

        
#'''!!!!!!****** start editing here to add module motif limit to two ******!!!!!!'''   
        
    while motifIndex < len(listOfSortedMotifs):
        combosInWindow = []
        
        # add motif that starts window
        combosInWindow.append(listOfSortedMotifs[motifIndex][0])
        
        # check to make sure the next motif is in the window and that it does not overlap the first motif
        #   if so: add combo to 
        motifCounter = 1
        
        while motifIndex + motifCounter < len(listOfSortedMotifs)\
        and listOfSortedMotifs[motifIndex+motifCounter][1] - listOfSortedMotifs[motifIndex][1] < windowLen \
        and not listOfSortedMotifs[motifIndex+motifCounter][1] - len(listOfSortedMotifs[motifIndex+motifCounter][0]):
            listIndex = motifIndex+motifCounter
            combosInWindow.append(listOfSortedMotifs[motifIndex+motifCounter][0])
            motifCounter += 1
        
        listOfMotifSets.append(Set(combosInWindow))   
        motifIndex += 1
    
    # use Set() to get nr set of motif sets then iterate through the sets generating all posible combinations 
    # of all lengths in each set(since the motif sets have already passed the windowLen test, all combinations generated from these sets
    # must also be no more than windowLen bp apart) 
    setOfMotifSets = Set(listOfMotifSets)
    
    for motifSet in setOfMotifSets:
        #  list-ify to prevent index shifting
        listFromSet = list(motifSet)
        
        numOfMotifs = len(listFromSet)
        
        #  to catch all combintiations of all lengths less than defined in loop for the seed set of motifs
        allCombinations = []
        comboLen = 1
        
        #  Loop to generate all combos of all lengths less than 4
        while comboLen <= numOfMotifs and comboLen <= 4:
            comboObj = Combination(listFromSet,comboLen)
            for item in comboObj:
                allCombinations.append(Set(item))
            comboLen+=1
        
        #  add list of generated combos to all sets found for this AGAP, we will remove redundancy later on
        listOfMotifSets += allCombinations
        
        
    #  Here by set-ifying listOfMotifSets we remove any redundant sets of motif
    ##leng = len(listOfMotifSets)    
    listOfMotifSets = Set(listOfMotifSets)
    ##leng = len(listOfMotifSets) 
    #  assign listOfMotifSets to its AGAP key as part of a list with index[0](listOfMotifSets) and index[1] representing 
    #  the presense<1> or absence<1> of a module to be searched in a later function
    dictOfMotifSetsByAGAP[AGAPname]= [list(listOfMotifSets), 0]
    w = 0
    
    #for each in listOfMotifSets:
        #print each
        
        
    print 'buildModulesForAGAP found %s modules for %s.' % (str(len(listOfMotifSets)), AGAPname)