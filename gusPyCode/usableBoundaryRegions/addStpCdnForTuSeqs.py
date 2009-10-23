

inFile  = '/Users/biggus/Documents/James/Data/Tu_miRNA/Cq_500afterCoding.newCoords.usuable.txt'
outFile = '/Users/biggus/Documents/James/Data/Tu_miRNA/Cq_500afterCoding.newCoords.usuable.stpCdn.txt'

inFile = map(lambda line: line.strip().split('\t'), open(inFile).readlines())

for i in range(0, len(inFile)):
    print inFile[i]
    if inFile[i][2] == '-1':
        inFile[i][4] = str(int(inFile[i][4])+3)
        inFile[i][5] = str(int(inFile[i][4])-int(inFile[i][3])+1)
    elif inFile[i][2] == '1':
        inFile[i][3] = str(int(inFile[i][3])-3)
        inFile[i][5] = str(int(inFile[i][4])-int(inFile[i][3])+1)

outFile = open(outFile,'w')
for x in inFile:
    outFile.write('\t'.join(x)+'\n')
    
    
print 'Done.'