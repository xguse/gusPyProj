from TAMO.seq import Fasta
from JamesDefs import revComp

fFile = '/Users/biggus/Documents/James/Data/Tu_miRNA/Fastas/Ag_500afterCoding.usuable.stpCdn.fas'
sFile = '/Users/biggus/Documents/James/Data/Tu_miRNA/miRNAs/testSeeds.fas'
oFile = '/Users/biggus/Documents/James/Data/Tu_miRNA/SeedCountOutPut/counts/testSeedsOn_Ag_500afterCoding.txt'

# --------- Fasta Prep ---------
fastas    = Fasta.file2dict(fFile)
seqNames  = fastas.keys()
seqNames.sort()
# seqs are softMasked.  This unMaskes them.
for seq in fastas:
    fastas[seq] = fastas[seq].upper()

# --------- Seed Prep ---------
seeds     = Fasta.file2dict(sFile)
seedNames = seeds.keys()
seedNames.sort()
# to make sure we are only looking for uppercase strings
for seed in seeds:
    seeds[seed] = seeds[seed].upper()


results = ['#seqName\t'+'\t'.join(seedNames)]

def findSeedsInSeq(seeds,seedNames,seqStr,seqName):
    '''take dict of seeds, a seq, and its name. Return tsv string of 
    seqName followed by 0s and 1s corelating with presence
    or absence of seed in seq'''
    
    hits = []
    for s in seedNames:
        if revComp(seeds[s]) in seqStr:
            hits.append('1')
        else:
            hits.append('0')
            
    return '%s\t%s' % (seqName,'\t'.join(hits))


print 'counting...'
for seq in seqNames:
    seqStr = fastas[seq]
    results.append(findSeedsInSeq(seeds,seedNames,seqStr,seq))

print 'writing to file...'    
oFile = open(oFile, 'w')
for l in results:
    oFile.write(l+'\n')
    
print 'Done.'