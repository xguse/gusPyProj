from TAMO.seq import Fasta
from gusPyCode.defs.bioDefs import geneList2FastaDict
from gusPyCode.defs.mosqData import promoterSeqPaths


geneList = map(lambda l: l.strip(), open('/Users/biggus/Documents/James/Collaborations/Campbell/data/CCupAt4Days.genes.txt', 'rU'))

sourceFasta = promoterSeqPaths.Aa_2000bpUp_hardMasked_shuf3

oFile = '/Users/biggus/Documents/James/Collaborations/Campbell/data/CCupAt4Days.masked.shuffled.3.fas'





newFasta = geneList2FastaDict(geneList, sourceFasta, hardMasked=True)

newFasta = Fasta.text(newFasta)
    
oFile = open(oFile, 'w')
oFile.write(newFasta)

print 'Done'