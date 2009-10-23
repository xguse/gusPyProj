from gusPyCode.defs.doug_hypergeometric import hyperGeoPvalue
from gusPyCode.defs import JamesDefs


#========================= User Defined Variables =========================
inFile  = '/Users/biggus/Documents/James/Data/OrthologDefs/nrOrthos/Data_Attempt2/MOW_Dm_Orthos_NonCompressed.txt'
outFile = '/Users/biggus/Documents/James/Data/OrthologDefs/nrOrthos/Data_Attempt2/MOW_Dm_Orthos.Compressed.txt'

currentGuideColumn = 3

#==========================================================================

#--------- Script Specific Function Definitions ---------------------

def extractAndCompress(unCompressedList, guideColumn):
    i = 0
    tempList = []    
    tempList.append(unCompressedList[i])
    i+=1
    
    # Copy matching records to tempList
    w = 'y'
    while w == 'y' and i < len(unCompressedList):
        if unCompressedList[i][guideColumn] == tempList[0][guideColumn]:
            tempList.append(unCompressedList[i])
            i+=1
        else:
            w = 'n'
    # Delete copied records so that next iteration starts at next set ALSO del i so it can be used again clean       
    del unCompressedList[0:len(tempList)], i
    
    # Compress the list of lists into one list with no blank indexes
    resultList = tempList[0]
    for i in range(0,len(tempList)):
        for j in range(0,len(tempList[i])):
            if tempList[i][j] != '':
                resultList[j] = tempList[i][j]
                
    # Convert resultList into a string
    resultString = resultList.pop(0)
    while resultList:
        resultString = '%s\t%s' % (resultString, resultList.pop(0))
    resultString+='\n'
    x=1
    return resultString
    
    


#--------------------------------------------------------------------






inFile = map(lambda line : line.rstrip('\n').split('\t'), open(inFile, 'rU').readlines())

compressedList = []

while inFile:
    compressedList.append(extractAndCompress(inFile,currentGuideColumn))

print "final length of list is: %i" % (len(compressedList))

outFile = open(outFile, 'w')
outFile.writelines(compressedList)

        