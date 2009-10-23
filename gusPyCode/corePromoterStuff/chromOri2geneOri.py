from bioDefs import revComp

# To convert Aa putative TSS seqs into +1 Ori based on gene's Ori


chromOriSeqs = map(lambda line: line.strip(), open('/Users/biggus/Documents/James/Data/AedesCorePromoterWork/output/Inr_DPE/Aa_Ensbl49_AaegL1.1.plus100minus100._InrDPE_.newInr.chromStrand.fa','rU').readlines())
outFile = '/Users/biggus/Documents/James/Data/AedesCorePromoterWork/output/Inr_DPE/Aa_Ensbl49_AaegL1.1.plus100minus100._InrDPE_.newInr.geneStrand.local.fa'
geneOriSeqs  = []


chromOriSeqs = zipped = zip(chromOriSeqs[:-1:2], chromOriSeqs[1::2])

for item in chromOriSeqs:
    fieldsList = item[0].split(' ')
    seq = item[1]
    if fieldsList[-1] == '(1:-1)':
        fieldsList[-1] = '(-1:1)'
        seq = revComp(seq)
        geneOriSeqs.append(' '.join(fieldsList)+'\n')
        geneOriSeqs.append(seq+'\n')
    else:
        geneOriSeqs.append(item[0]+'\n')
        geneOriSeqs.append(item[1]+'\n')
        
        
outFile = open(outFile, 'w')
outFile.writelines(geneOriSeqs)
    
print 'Done.'

