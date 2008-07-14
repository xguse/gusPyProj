import sys
import time
import JamesDefs

#--------- Script Specific Function Definitions ---------------------
def resolver(groupedList, resolverArgs):
    
    returnList = []
    tab = '\t'
    for targetGeneConflicts in groupedList:
        # due to a current cluster kludge of a solution we must check for a "kludge group"
        # which should be the final rec
        if targetGeneConflicts[0][0] == 'kludge':
            break # should be the last group anyway
        
        # open clean lists for string and number versions of coords
        strCoordsList = []
        intCoordsList = []    
        for conflictRec in targetGeneConflicts:
            # int('') throws a TypeError so we must filter out the coords to a clean list to convert them from str to int
            if conflictRec[5] != '':
                strCoordsList.append(conflictRec[resolverArgs['conflictRegionStrt']])
            if conflictRec[6] != '':
                strCoordsList.append(conflictRec[resolverArgs['conflictRegionEnd']])
                
        # Convert str coords into int coords
        for strCoord in strCoordsList:
            intCoordsList.append(int(strCoord))
            
        # decide which scheme to use based on Up or Dwn stream coords
        if resolverArgs['whichBoundary'] == 'upStream':
            strandList = [['+','1'],['-','-1']]
        elif resolverArgs['whichBoundary'] == 'downStream':
            strandList = [['-','-1'],['+','1']]
            
        if targetGeneConflicts[0][resolverArgs['strandField']]  in strandList[0]:
            # calculate usable portion of conflicted gene boundaries and add record to the returnList
            
                                        
            lowerBoundProximal     = int(targetGeneConflicts[0][resolverArgs['lowerBoundProximal']])
            maxConflictCoordPlus1  = max(intCoordsList)+1
            useableRegLen          = lowerBoundProximal - maxConflictCoordPlus1 + 1
            
            #                  geneID                                chromo                                     strand                          
            usableBdryReg = targetGeneConflicts[0][1] +tab+ targetGeneConflicts[0][2] +tab+ targetGeneConflicts[0][resolverArgs['strandField']] +tab+ str(maxConflictCoordPlus1) +tab+ str(lowerBoundProximal) +tab+ str(useableRegLen) +'\n'
            returnList.append(usableBdryReg)
            usableBdryReg = None  # clear variable for next time
        
        elif targetGeneConflicts[0][resolverArgs['strandField']]  in strandList[1]:
            # calculate usable portion of conflicted gene boundaries and add record to the returnList
            
                                        
            higherBoundProximal    = int(targetGeneConflicts[0][resolverArgs['higherBoundProximal']])
            minConflictCoordMinus1 = min(intCoordsList)-1
            useableRegLen          = minConflictCoordMinus1 - higherBoundProximal + 1
            
            #                  geneID                               chromo                                     strand                          
            usableBdryReg = targetGeneConflicts[0][1] +tab+ targetGeneConflicts[0][2] +tab+  targetGeneConflicts[0][resolverArgs['strandField']] +tab+ str(higherBoundProximal) +tab+ str(minConflictCoordMinus1) +tab+ str(useableRegLen) +'\n'
            returnList.append(usableBdryReg)
            usableBdryReg = None  # clear variable for next time
        else:
            sys.exit('StrandField has unexpected value for record '+targetGeneConflicts[0][1])
    
    return returnList
#--------------------------------------------------------------------


############### list of user variables ####################

# file to import from
inFile  = '/Users/biggus/Documents/James/Data/2KB/2kb_Sequence/2kb_Drosophila/2KBdown_drosophilaBoundaryCoordsTSS.both.fjoin.txt'

# file to write output to
outFile = '/Users/biggus/Documents/James/Data/2KB/2kb_Sequence/2kb_Drosophila/2KBdown_drosophilaBoundaryCoordsTSS.both.fjoinResolved.txt '

whichBoundary = 'downStream' # 'upStream' or 'downStream'




# open and create handle for inFile
conflictFile = open(inFile, 'r')

# open and create handle for outFile
resFile = open(outFile, 'w')
tick = time.time()

# read file into list 
conflictList = conflictFile.readlines()
# remove trailing '\n' from every record
LEN_cL = len(conflictList)
i = 0
while i < LEN_cL:
    conflictList[i] = conflictList[i].rstrip('\n')
    i = i + 1


# group file by target gene id using groupByField 
fjoinOutByGeneIDList = JamesDefs.groupByField(conflictList, 1)

resolverArgs = {
                    'strandField' : 4,
                    'lowerBoundProximal' : 10,
                    'higherBoundProximal' : 11,
                    'conflictRegionStrt' : 18,
                    'conflictRegionEnd' : 19,
                    'whichBoundary':whichBoundary
                }

resolvedBoundariesList = resolver(fjoinOutByGeneIDList, resolverArgs)

resFile.writelines(resolvedBoundariesList)

tock = time.time()
print 'The process took ', tock - tick, ' seconds.'