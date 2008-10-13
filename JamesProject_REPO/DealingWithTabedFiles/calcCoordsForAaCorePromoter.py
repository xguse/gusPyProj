

inFile  = '/Users/biggus/Documents/James/Data/AedesCorePromoterWork/Aa_Ensbl49_AaegL1.1.NR.txt'
outFile = '/Users/biggus/Documents/James/Data/AedesCorePromoterWork/Aa_Ensbl49_AaegL1.1.PlusMinus50coords.txt'


inFile = map(lambda line: line.strip().split('\t'), open(inFile, 'rU').readlines())


passedRecs = 0

def calcCoords(line, passedRecs):
    if int(line[6]) == 1:
        if int(line[7])-50 < 0:
            passedRecs+=1
            print "%s was too short." % (line[0])
        else:            
            return '%s\t%s\t%s\t%s\t%s\t%s\n' % (line[3],line[0],line[1],line[6],int(line[7])-50,int(line[7])+50)
    elif int(line[6]) == -1:
        if int(line[8])-50 < 0:
            passedRecs+=1
            print "%s was too short." % (line[0])
        else:
            return '%s\t%s\t%s\t%s\t%s\t%s\n' % (line[3],line[0],line[1],line[6],int(line[8])-50,int(line[8])+50)

coordList = []
for line in inFile:
    coordList.append(calcCoords(line,passedRecs))
    
outFile = open(outFile,'w')

outFile.writelines(coordList)
