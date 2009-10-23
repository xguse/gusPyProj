coordFile  = map(lambda line: line.strip().split('\t'), open('/Users/biggus/Documents/James/Data/AedesCorePromoterWork/sourceCoords/Aa_Ensbl49_AaegL1.1.NR.txt','rU').readlines())
queryFile  = map(lambda line: line.strip().split('\t'), open('/Users/biggus/Documents/James/Data/AedesCorePromoterWork/output/TATA_Inr/Aa_Ensbl49_AaegL1.1.plus50minus100._TataInr_.InrCoords.txt','rU').readlines())
outFile    = '/Users/biggus/Documents/James/Data/AedesCorePromoterWork/output/TATA_Inr/Aa_Ensbl49_AaegL1.1.plus50minus100._TataInr_.newInrCoords.txt'
newInrFile = []


for qLine in queryFile:
    for cLine in coordFile:
        if qLine[0] == cLine[3]:
            # For Plus Strand genes
            if cLine[6] == '1':
                name     = cLine[3]
                chrom    = cLine[1]
                strand   = cLine[6]
                newInr   = int(cLine[7])-(100-int(qLine[1])) # 100 is where old slice are centered so this corrects qLine[1]
                newStart = newInr - 100
                newStop  = newInr + 100
                newInrFile.append('%s\t%s\t%s\t%s\t%s\t%s\n' % (name,chrom,strand,newInr,newStart,newStop))
                ##x=1
            if cLine[6] == '-1':
                name     = cLine[3]
                chrom    = cLine[1]
                strand   = cLine[6]
                newInr   = int(cLine[8])+(100-int(qLine[1])) # 100 is where old slice are centered so this corrects qLine[1]
                newStart = newInr - 100
                newStop  = newInr + 100
                newInrFile.append('%s\t%s\t%s\t%s\t%s\t%s\n' % (name,chrom,strand,newInr,newStart,newStop))
                ##if cLine[3] == 'AAEL014559-RA':
                    ##x=1

print '%s adjustments made.' % (len(newInrFile))

outFile = open(outFile, 'w')

outFile.writelines(newInrFile)
print 'Done.'


