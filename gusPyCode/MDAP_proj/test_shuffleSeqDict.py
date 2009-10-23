from gusPyCode.MDAP_proj.MDAP_defs import shuffleSeqDict
from TAMO.seq import Fasta
from gusPyCode.defs.bioDefs import softMaskDict2HardMask
from time import time
from gusPyCode.defs.mosqData import promoterSeqPaths

# User Variables:
inFile   = promoterSeqPaths.Aa_2000bpUp_softMasked
outFile  = '/Users/biggus/Documents/James/Data/2KB/2kb_Sequence/2kb_Aedes/aedes2KBupStreamTSS.UnMasked.geneStrand.shuffledSeqs.1.fas'
hardMask = None


d = Fasta.load(inFile)
#d = {1:'AACTGCANACTGACNNNACTGATGNNN'}

if not hardMask:
    for x in d:
        d[x] = d[x].upper()

t1 = time()
sD = shuffleSeqDict(d)
t2 = time()

Fasta.write(sD,outFile)

print 'Shuffling took %.2f min.' % ((float(t2)-t1)/60)