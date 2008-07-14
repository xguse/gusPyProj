import sys
import time
############### list of user variables ####################

# file to import from
inFile = '/Users/biggus/Documents/MBGB/Rotations/James/Data/Sequence/Anopheles/testCodingBoundsOut/anophCodingCoordsProtCod_fjoin_TEST.txt'

# file to write output to
outFile = '/Users/biggus/Documents/MBGB/Rotations/James/Data/Sequence/Anopheles/testCodingBoundsOut/anophCodingCoordsProtCod_fjoin_TEST_RESOLVED.txt'

strandField = 4

lowerBoundProximal = 10
higherBoundProximal = 11

conflictRegionStrt = 18
conflictRegionEnd = 19


# open and create handle for inFile
conflictFile = open(inFile, 'r')

# open and create handle for outFile
resFile = open(outFile, 'w')


# -read in data line by line testing for a change in upStrm name(col:1 of line:0 ) in fjoin outFile

# -create list of lists with 1-ary list housing lines and 2-ary lists holding parsed fields


resolvedCoords = []    # to catch resolved coordinates

upStreamCons = []    # working home of lines representing upStream conflicts

tick = time.clock()

upStreamCons.append(conflictFile.readline().split('\t')) # seed with first line of file

nextLineFields = conflictFile.readline().split('\t')

while nextLineFields[0] != '':

    if nextLineFields[1] == upStreamCons[0][1]: # if upStream name matches seeded line 
        upStreamCons.append(nextLineFields)    #then add the split line to list
        nextLineFields = conflictFile.readline().split('\t') # and advance to next line of file

    else:  # if name doesnt match, then start processing 'upStreamCons' list
        if upStreamCons[0][strandField] in ['+','1']:
            # initiate clean vars each round to prevent 'old' data hanging around
            #print '1 Pos\n'
            transStrt =''
            conflictingCoords = []
            proxCfltEdge =''
            upStrmLen = 0
            transStrt = int(upStreamCons[0][lowerBoundProximal])    # 
            for each in upStreamCons:    # append all start/stop coords of exons that conflict with upStrm 1K
                #print '2 Pos\n'
                conflictingCoords.append(int(each[conflictRegionStrt]))
                conflictingCoords.append(int(each[conflictRegionEnd]))
            proxCfltEdge = max(conflictingCoords)    # most proximal to translationStartSite is MAX for PLUS strand
            upStrmLen = transStrt - proxCfltEdge + 1    # compute length of upStrm Feature
            resolvedCoords.append(upStreamCons[0][1]+'\t'+upStreamCons[0][strandField]+'\t'+str(min(transStrt,proxCfltEdge))+'\t'+str(max(transStrt,proxCfltEdge))+'\t'+str(upStrmLen))
            upStreamCons = []    # clear var
            upStreamCons.append(nextLineFields)    # seed with current line of file
        elif upStreamCons[0][strandField] in ['-','-1']:
            # initiate clean vars each round to prevent 'old' data hanging around
            #print '1 Neg\n'
            transStrt =''
            conflictingCoords = []
            proxCfltEdge =''
            upStrmLen = 0
            transStrt = int(upStreamCons[0][higherBoundProximal])    # 
            for each in upStreamCons:    # append all start/stop coords of exons that conflict with upStrm 1K
                #print '2 Neg\n'
                conflictingCoords.append(int(each[conflictRegionStrt]))
                conflictingCoords.append(int(each[conflictRegionEnd]))
            proxCfltEdge = min(conflictingCoords)    # most proximal to translationStartSite is MIN for MINUS strand
            upStrmLen = proxCfltEdge - transStrt + 1    # compute length of upStrm Feature
            resolvedCoords.append(upStreamCons[0][1]+'\t'+upStreamCons[0][strandField]+'\t'+str(min(transStrt,proxCfltEdge))+'\t'+str(max(transStrt,proxCfltEdge))+'\t'+str(upStrmLen))
            upStreamCons = []    # clear var
            upStreamCons.append(nextLineFields)    # seed with current line of file
        else:
            sys.exit('I did not see +/1 OR -/-1 when checking for strand')


resFile.write('\n'.join(resolvedCoords))

tock = time.clock()

conflictFile.close()
resFile.close()

print 'The process took ', tock - tick, ' seconds.'