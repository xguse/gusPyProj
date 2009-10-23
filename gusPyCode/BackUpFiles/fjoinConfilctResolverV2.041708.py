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
                strCoordsList.append(conflictRec[18])
            if conflictRec[6] != '':
                strCoordsList.append(conflictRec[19])
                
        # Convert str coords into int coords
        for strCoord in strCoordsList:
            intCoordsList.append(int(strCoord))
            
        if targetGeneConflicts[0][resolverArgs['strandField']]  in ['+','1']:
            # calculate usable portion of conflicted gene boundaries and add record to the returnList
            
                                        
            useableRegStrt = min(int(targetGeneConflicts[0][resolverArgs['lowerBoundProximal']]),max(intCoordsList)+1)
            useableRegEnd  = max(int(targetGeneConflicts[0][resolverArgs['lowerBoundProximal']]),max(intCoordsList)+1)
            useableRegLen  = useableRegEnd - useableRegStrt + 1
            
            #                  geneID                                        strand                          
            usableBdryReg = targetGeneConflicts[0][1] +tab+ targetGeneConflicts[0][resolverArgs['strandField']] +tab+ str(useableRegStrt) +tab+ str(useableRegEnd) +tab+ str(useableRegLen) +tab+ '\n'
            returnList.append(usableBdryReg)
            usableBdryReg = None  # clear variable for next time
        
        elif targetGeneConflicts[0][resolverArgs['strandField']]  in ['-','-1']:
            # calculate usable portion of conflicted gene boundaries and add record to the returnList
            
                                        
            useableRegStrt = min(int(targetGeneConflicts[0][resolverArgs['higherBoundProximal']]),min(intCoordsList)-1)
            useableRegEnd  = max(int(targetGeneConflicts[0][resolverArgs['higherBoundProximal']]),min(intCoordsList)-1)
            useableRegLen  = useableRegEnd - useableRegStrt + 1
            
            #                  geneID                                         strand                          
            usableBdryReg = targetGeneConflicts[0][1] +tab+ targetGeneConflicts[0][resolverArgs['strandField']] +tab+ str(useableRegStrt) +tab+ str(useableRegEnd) +tab+ str(useableRegLen) +tab+ '\n'
            returnList.append(usableBdryReg)
            usableBdryReg = None  # clear variable for next time
        else:
            sys.exit('StandField has unexpected value for record '+targetGeneConflicts[0][1])
    
    return returnList

############### list of user variables ####################

# file to import from
inFile = '/Users/biggus/Documents/MBGB/Rotations/James/Data/Sequence/Anopheles/TestingMyPipe/anophCodingCoords_fjoin_BOTH.out.sorted.txt'

# file to write output to
outFile = '/Users/biggus/Documents/MBGB/Rotations/James/Data/Sequence/Anopheles/TestingMyPipe/anophCodingCoords_fjoin_BOTH.resolved.txt'





# open and create handle for inFile
conflictFile = open(inFile, 'r')

# open and create handle for outFile
resFile = open(outFile, 'w')
tick = time.clock()

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
                    'whichBoundary':'upStream'
                }

resolvedBoundariesList = resolver(fjoinOutByGeneIDList, resolverArgs)

resFile.writelines(resolvedBoundariesList)

tock = time.clock()
print 'The process took ', tock - tick, ' seconds.'