import JamesDefs
## better version in defs_moduleByTableLoopUp
#def makeFwdAndRevCompRegExObj(fwdAndRevCompRegExStrings):
    #import re
    #fwdAndRevCompRegExStrings = '(%s|%s)' % (fwdAndRevCompRegExStrings[0], fwdAndRevCompRegExStrings[1])
    
    #fwdRevComp_regExObj = re.compile(fwdAndRevCompRegExStrings, re.IGNORECASE)
    #return fwdRevComp_regExObj

def createSearchModules(listOfMotifs,mostMotifsInModule):
    from time import time
    from probstat import Combination
    
    
    t1= time()
    listOfSearchModules = []
    len_listOfMotifs = len(listOfMotifs)
    
    comboLen = 1
    while comboLen <= len_listOfMotifs and comboLen <= mostMotifsInModule:
        comboObj = Combination(listOfMotifs,comboLen)
        for item in comboObj:
            listOfSearchModules.append(item)
        comboLen+=1
    
    len_listOfSearchModules = len(listOfSearchModules)
    t2 = time()
    print "Returning %s search modules in %s sec." % (str(len_listOfSearchModules), str(t2-t1))
    return listOfSearchModules

#==============

def countModuleInAll(module, seqDict):
    
    
    
    import re
    
    ##motifList = moduleStrs.split('\t')
    
    #  Calculate revComp for each motif, convert them into perl-like regular expressions, compile those RegExStrings into
    #    python regEx objs and modify data structure to be list of lists with each 
    #    2ary list = [regExObj_Fwd, regExObj_Rev]
   
    i = 0
    for IUPACmotif in module:
        #  convert data struct and calc revComp
        module[i] = [IUPACmotif, [IUPACmotif, JamesDefs.revComp(IUPACmotif)]]
        c = 0
        for each in module[i][1]:
            module[i][1][c] = JamesDefs.iupac2regex(module[i][1][c])
            c+=1
        module[i][1] = makeFwdAndRevCompRegExObj(module[i][1])
        i+=1

    
    
    
    
    
    
    
    
    #  Loop over list and count module in fwd and revComp oris for presence absense
    #  Sum total hits in totalHits
    
    totalHits = 0
    
    for record in seqDict:
        #  initiate modulePresent attrib to 0
        seqDict[record].modulePresent = 0
                
        hit = findModuleInLength(module, seqDict[record], 500)
        totalHits+=hit
        
    return totalHits

def convertMotifList(motifList):
    i = 0
    while i < len(motifList):
        motifList[i] = [motifList[i], JamesDefs.iupac2regex(motifList[i])]
        i += 1

def findModuleInLength(listOfModuleInstances, SeqObj, windowLen):
    seqToSearch = SeqObj.seq.tostring()
    
    
    i = 0    
    while i < len(seqToSearch):
        #  pull out window to search
        currentWindow = seqToSearch[i:i+windowLen]
        
        # search for motifs 
        searchResult = findModuleInWindow(listOfModuleInstances, currentWindow, 0)
        
        if searchResult == 'success':
            SeqObj.modulePresent = 1
            
            #  If search succeded stop search and move to next module
            return 1
        else:
            i+=1
            
    return 0 

def findModuleInWindow(listOfModuleInstances, windowSeq, moduleMember=0): # module member MUST start at 0
    """
      Search windowSeq with fwd regex object
      If successful, call next instance of findModuleInWindow() with next module member
      If unsuccessful with Fwd, try RevComp
      If successful call next instance of findModuleInWindow() with next module member
      If unsuccessful return 'failed'
    """
    #  Extract compiled regEx versions of each 
    
    #  Test for module member
    if listOfModuleInstances[moduleMember][1].search(windowSeq):
        # -check to see if this is last member
        if moduleMember+1 == len(listOfModuleInstances):
            # --if last memeber return success
            return 'success'
        else:
            # --if not last member call next instance of findModuleInWindow() with next module member 
            moduleMember+=1
            result = findModuleInWindow(listOfModuleInstances, windowSeq, moduleMember)
            return result
        
    else:
        return 'failed'
            