print 'loading seq file...'
seqFile    = map(lambda line: line.strip(), open('/Users/biggus/Documents/James/Data/AedesCorePromoterWork/output/TATA_Inr/Aa_Ensbl49_AaegL1.1.plus100minus100._TataInr_.newInr.geneStrand.fa','rU').readlines())
outFile    = '/Users/biggus/Documents/James/Data/AedesCorePromoterWork/output/TATA_Inr/Aa_Ensbl49_AaegL1.1.plus100minus100._TataInr_.newInr.geneStrand.nrSeq.fa'
fastaList = zipped = zip(seqFile[:-1:2], seqFile[1::2])

nrSeqs   = []
nrFastas = []

namesOfExcludedSeqs = []

for fasta in fastaList:
    if fasta[1] not in nrSeqs:
        nrSeqs.append(fasta[1])
        nrFastas.append(fasta)
    else:
        namesOfExcludedSeqs.append(fasta[0])
        
        
#open outFile and write data/close file
outFile = open(outFile, 'w')

for each in nrFastas:
    outFile.write('>%s\n%s\n'%(each[0],each[1]))



print """%s out of %s fasta entries included in nr file.
The following are headers for those seqs excluded:\n""" % (len(nrFastas),len(fastaList))
for name in namesOfExcludedSeqs:
    print name