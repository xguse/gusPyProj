############### list of user variables ####################

# file to import from
rFile = '/Users/biggus/Documents/James/Data/OrthologDefs/nrOrthos/Data_Attempt2/test4Mosq3WayOrthos.nr.txt'

# file to write output to
outFile = '/Users/biggus/Documents/James/Data/OrthologDefs/nrOrthos/Data_Attempt2/test4Mosq3WayOrthos.nrFirstIndex.txt'

# open and read data into var
rFile = map(lambda line : line.rstrip('\n').split('\t'), open(rFile, 'rU').readlines())



# invoke list to hold results
nrList = []
firstIndexes = []
secondIndexes = []
thirdIndexes = []
repeatList = []

# read one line from file at a time and:
# -strip the trailing \n
# -ask if it is in the results list
# -if no put it there, if yes go to next line

for line in rFile:
    
    if line[0] not in firstIndexes and line[1] not in secondIndexes and line[2] not in thirdIndexes:
        nrList.append("\t".join(line))
        firstIndexes.append(line[0])
        secondIndexes.append(line[1])
        thirdIndexes.append(line[2])
   
    else:
        if line[0] not in repeatList:
            repeatList.append("\t".join(line))
        


#open outFile and write data/close file
nrFile = open(outFile, 'w')
nrFile.write('\n'.join(nrList))

print str(len(nrList))+' non-repetative sequences written. '+str(len(rFile)-len(nrList))+' items excluded.'
print 'The following is a non-redundant list of those records that were repeated in original file:'
for everyOne in repeatList:
    print everyOne