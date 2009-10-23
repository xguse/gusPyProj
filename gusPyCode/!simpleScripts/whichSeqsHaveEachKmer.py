from TAMO.MotifTools import top_nmers,Motif
from TAMO import MotifTools
from TAMO.seq import Fasta
from gusPyCode.defs.bioDefs import revComp

seqFile     = '/Users/biggus/Documents/James/Collaborations/Campbell/data/mainTwoGenes.fas'
kmerFile    = '/Users/biggus/Documents/James/Collaborations/Campbell/data/mainTwoGenes.7mersInAll.tamoVers.txt'
#outFile     = '/Users/biggus/Documents/James/Collaborations/Campbell/data/mainTwoGenes.7mersInAll.tamoVers.txt'


seqs     = Fasta.file2dict(seqFile)
seqNames = sorted(seqs.keys())
kmers    = map(lambda l: l.strip(), open(kmerFile, 'rU').readlines())
kmers.append('AAHRRSSSSSSSSSMMMMM')
results = []
for kmer in kmers:
    print '>%s:' % (kmer)
    results.append('%s:' % (kmer))
    for seq in seqNames:
        if seqs[seq].find(kmer) != -1 or seqs[seq].find(revComp(kmer)) != -1:
            print '\t'+seq
            results.append(seq)
        else:
            print '---NOT FOUND---'
    
#outFile = open(outFile, 'w')
#for each in inAllSeqs:
    #outFile.write(each+'\n')




print 'Done'