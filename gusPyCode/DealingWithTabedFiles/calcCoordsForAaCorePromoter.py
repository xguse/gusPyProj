

inFile  = '/Users/biggus/Documents/James/Data/AedesCorePromoterWork/sourceCoords/Aa_Ensbl49_AaegL1.1.NR.txt'
outFile = '/Users/biggus/Documents/James/Data/AedesCorePromoterWork/sourceCoords/Aa_Ensbl49_AaegL1.1.Plus50Minus100coords.txt'


inFile = map(lambda line: line.strip().split('\t'), open(inFile, 'rU').readlines())


passedRecs = 0

def calcCoords(line, passedRecs):
    if int(line[6]) == 1:
        return '%s\t%s\t%s\t%s\t%s\t%s' % (line[3],line[0],line[1],line[6],int(line[7])-100,int(line[7])+100)
    elif int(line[6]) == -1:
        return '%s\t%s\t%s\t%s\t%s\t%s' % (line[3],line[0],line[1],line[6],int(line[8])-100,int(line[8])+100)

coordList = []
for line in inFile:
    coordList.append(calcCoords(line,passedRecs))

# check for neg numbers in the coords and add trailing \n
c = 0
for i in range(0, len(coordList)):
    coordList[i] = coordList[i].split('\t')
    if int(coordList[i][4]) < 0:
        print '%s was too short.' % (coordList[i][0])
        coordList[i] = ''
        c+=1
    else:
        coordList[i] = "\t".join(coordList[i])+'\n'

outFile = open(outFile,'w')
for each in coordList:
    #print each
    outFile.write(each)

#outFile.writelines(coordList)

print '--------\nDone.\n%s seqs were too short and were ommited.\n**see above**' % (c)