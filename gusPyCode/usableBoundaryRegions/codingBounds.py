import sys
import time
from gusPyCode.defs import JamesDefs


#--------- Script Specific Function Definitions ---------------------
def combineExons(groupedList, BdryLen):
    
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
            if strCoord.startswith('#'):
                continue
            intCoordsList.append(int(strCoord))
        
        # Write new one line record for the gene
        
        # Correct any lowBdryLenStrt/lowBdryLenEnd coords that are negative 
        #
        # Start with lowBdryLenStrt:
        correctedLowBdryLenStrt = None # initialize clean var
        if min(intCoordsList)-BdryLen < 1:
            correctedLowBdryLenStrt = 1
        else:
            correctedLowBdryLenStrt = min(intCoordsList)-BdryLen
        
        # Now lowBdryLenEnd    
        correctedLowBdryLenEnd = None # initialize clean var
        if min(intCoordsList)-1 < 1:
            correctedLowBdryLenEnd = 1
        else:
            correctedLowBdryLenEnd = min(intCoordsList)-1

        # Format outPut record    
        #                AGAP#_0        Chromo_1        bioType_2       Strand_3            geneStrt_4                     geneEnd_5                 NumOfExons_6        
        oneLineRecord = gene[0][0]+'\t'+gene[0][1]+'\t'+gene[0][2]+'\t'+gene[0][3]+'\t'+str(min(intCoordsList))+'\t'+str(max(intCoordsList))+'\t'+str(len(gene))+'\t'+ \
                     str(max(intCoordsList)-min(intCoordsList)+1)+'\t'+str(correctedLowBdryLenStrt)+'\t'+str(correctedLowBdryLenEnd)+'\t'+str(max(intCoordsList)+1)+'\t'+str(max(intCoordsList)+BdryLen)+'\n'
                   #           ChromoCoverage_7                            lowBdryLenStrt_8                   lowBdryLenEnd_9                   hiBdryLenStrt_10                   hiBdryLenEnd_11

        returnList.append(oneLineRecord)
    return returnList
#--------------------------------------------------------------------

#========================= User Defined Variables =========================

# File paths:
sourceFile = '/Users/biggus/Documents/James/Data/Tu_miRNA/MegyDBdumpData/Cq_CpipJ1-2_GeneTranscrExon_112408.noBlanks.txt'
boundaryOutPutFile = '/Users/biggus/Documents/James/Data/Tu_miRNA/Cq_500afterCoding.newCoords.txt'
    
#  Boundary Region Length (1000 or 3000 bp etc)
bdryLen = 500

#==========================================================================



srcFile   = open(sourceFile, 'r')
boundaryFile  = open(boundaryOutPutFile, 'w')




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
oneLineRecordList = combineExons(groupedList, bdryLen)

# Write out oneLineRecordList to outFile
boundaryFile.writelines(oneLineRecordList)
boundaryFile.close()




print 'Tada!'







    