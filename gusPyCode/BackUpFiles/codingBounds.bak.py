import sys
import time
import JamesDefs
from fjoin import FJoin

#--------- Script Specific Function Definitions ---------------------
def combineExons(groupedList):
    
    returnList = []

    for gene in groupedList:
        # open clean lists for string and number versions of coords
        strCoordsList = []
        intCoordsList = []    
        for exon in gene:
            # int('') throws a TypeError so we must filter out the coords to a clean list to convert them from str to int
            if exon[5] != '':
                strCoordsList.append(exon[5])
            if exon[6] != '':
                strCoordsList.append(exon[6])
        
        # Convert str coords into int coords
        for strCoord in strCoordsList:
            intCoordsList.append(int(strCoord))
        
        # Write new one line record for the gene
        
        # Correct any low1000Strt/low1000End coords that are negative 
        #
        # Start with low1000Strt:
        correctedLow1000Strt = None # initialize clean var
        if min(intCoordsList)-1000 < 1:
            correctedLow1000Strt = 1
        else:
            correctedLow1000Strt = min(intCoordsList)-1000
        
        # Now low1000End    
        correctedLow1000End = None # initialize clean var
        if min(intCoordsList)-1 < 1:
            correctedLow1000End = 1
        else:
            correctedLow1000End = min(intCoordsList)-1

        # Format outPut record    
        #                AGAP#_0        Chromo_1        bioType_2       Strand_3            geneStrt_4                     geneEnd_5                 NumOfExons_6        
        oneLineRecord = gene[0][0]+'\t'+gene[0][1]+'\t'+gene[0][2]+'\t'+gene[0][3]+'\t'+str(min(intCoordsList))+'\t'+str(max(intCoordsList))+'\t'+str(len(gene))+'\t'+ \
                     str(max(intCoordsList)-min(intCoordsList)+1)+'\t'+str(correctedLow1000Strt)+'\t'+str(correctedLow1000End)+'\t'+str(max(intCoordsList)+1)+'\t'+str(max(intCoordsList)+1000)+'\n'
                   #           ChromoCoverage_7                            low1000Strt_8                   low1000End_9                   hi1000Strt_10                   hi1000End_11

        returnList.append(oneLineRecord)
    return returnList
#--------------------------------------------------------------------

#========================= User Defined Variables =========================

# File paths:
sourceFile = '/Users/biggus/Documents/MBGB/Rotations/James/Data/Sequence/Culex/Culex_Exon_Coding_Reordered.txt'
boundaryOutPutFile = '/Users/biggus/Documents/MBGB/Rotations/James/Data/Sequence/Culex/culexProcessed/culexCodingBoundsOut.txt'
#fjoinOutPutFile = '/Users/biggus/Documents/MBGB/Rotations/James/Data/Sequence/Anopheles/testCodingBoundsOut/anophCodingCoordsProtCod_fjoin_TEST.txt'
#finalOutPutFile = '/Users/biggus/Documents/MBGB/Rotations/James/Data/Sequence/Anopheles/testCodingBoundsOut/anophCodingCoordsProtCod_final_TEST.txt'

## fjoin Paths
#path2fjIn1 = boundaryOutPutFile
#path2fjIn2 = '/Users/biggus/Documents/MBGB/Rotations/James/Data/Sequence/Anopheles/TestingMyPipe/anophCodingCoords.bak.clean'

## other fjoin arguments: [will need to change 5 & 6 when doing ]
#columns1 = '2,5,6'  # (for input 1) Modify numbers to represent HUMAN column numbers for chromosome, start, and end, respectively
#columns2 = '2,6,7'  # (for input 1) See above
#==========================================================================



srcFile   = open(sourceFile, 'r')
boundaryFile  = open(boundaryOutPutFile, 'w')
#fjoinFile = open(fjoinOutPutFile, 'w')
#finalFile = open(finalOutPutFile,'w')




# Read source data into list
bioMartList = srcFile.readlines()

# remove trailing '\n' from every record
LEN_bML = len(bioMartList)
i = 0
while i < LEN_bML:
    bioMartList[i] = bioMartList[i].rstrip('\n')
    i = i + 1

# Grouping records by gene name and splitting record fields into lists
groupedList = JamesDefs.groupByField(bioMartList, 0)

# Combine exon records into a single gene line record with start and stop coords for coding region
# TranscriptID field will be removed and fields representing the number of exons encountered and 
# the chromosomal coverage will be appended respectivly to the end of each record
oneLineRecordList = combineExons(groupedList)

# Write out oneLineRecordList to outFile
boundaryFile.writelines(oneLineRecordList)
boundaryFile.close()

## Set up FJoin arguments and call FJoin
#fjoinArgs = ["-1", path2fjIn1, 
             #"-2", path2fjIn2, 
             #"-s", "both", 
             #"-o", fjoinOutPutFile,
             #"--columns1", columns1,
             #"--columns2", columns2]

#FJoin(fjoinArgs).go()



print 'Tada!'







    