##!/usr/bin/env python



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
        


#==================================

def dev_buildMotifPairsForAGAP(AGAPentryOfSortedMotifs, dictOfMotifSetsByAGAP, windowLen=500):
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
        #  list-ify to prevent index shifting
        listFromSet = list(motifSet)
        
        numOfMotifs = len(listFromSet)
        
        #  to catch all combintiations of 2
        combosOfTwo = []
        comboLen = 1
        
        #  Loop to generate all combos of all lengths less than 4
        if len(listFromSet) >= 2:
            comboObj = Combination(listFromSet,2)
            for item in comboObj:
                combosOfTwo.append(Set(item))
      
        
        #  add list of generated combos to all sets found for this AGAP, we will remove redundancy later on
        listOfMotifSets += combosOfTwo
        
        
    #  Here by set-ifying listOfMotifSets we remove any redundant sets of motif
    ##leng = len(listOfMotifSets)    
    listOfMotifSets = Set(listOfMotifSets)
    ##leng = len(listOfMotifSets) 
    #  assign listOfMotifSets to its AGAP key as part of a list with index[0](listOfMotifSets) and index[1] representing 
    #  the presense<1> or absence<0> of a module to be searched in a later function
    dictOfMotifSetsByAGAP[AGAPname]= [list(listOfMotifSets), 0]
    w = 0
    
    #for each in listOfMotifSets:
        #print each
        
        
    print 'buildModulesForAGAP found %s modules for %s.' % (str(len(listOfMotifSets)), AGAPname)























