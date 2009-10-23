from TAMO.MotifTools import top_nmers,Motif
from TAMO import MotifTools
from TAMO.seq import Fasta
from bioDefs import ifKmerInAll

seqFile     = '/Users/biggus/Documents/James/Collaborations/Campbell/data/mainTwoGenes.fas'
outFile     = '/Users/biggus/Documents/James/Collaborations/Campbell/data/mainTwoGenes.8mersInAll.txt'
kmerSize    = 8
scoreThresh = 0.999999

seqs = Fasta.file2dict(seqFile)



# create new dict to store the seqs' kmers
seqsKmers = {}
for i in seqs:
    seqsKmers[i] = top_nmers(kmerSize,[seqs[i]], purge_Ns = 1)   # for some reason top_nmers fails silently if given str instead of list

inAllSeqs = []
count = 0
for seq in seqsKmers:
    for kmer in seqsKmers[seq]:
        if ifKmerInAll(kmer,seqs,scoreThresh):
            if kmer not in inAllSeqs:
                inAllSeqs.append(kmer)
                count+=1
                print count


outFile = open(outFile, 'w')
for each in inAllSeqs:
    outFile.write(each+'\n')




print 'Done'