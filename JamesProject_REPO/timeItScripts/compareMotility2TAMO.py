import timeit
from time import time
from TAMO.MotifTools import Motif
from TAMO.seq import Fasta
import motility
from TAMO.MotifTools import load as loadTMOs

tOut = '/Users/biggus/Documents/James/Data/ReClustering/Python_CRM/tamoTimeIt.6memeMotifs.35seqs.30runs.txt'
mOut = '/Users/biggus/Documents/James/Data/ReClustering/Python_CRM/motilityTimeIt.6memeMotifs.35seqs.30runs.txt'

genes = 35
runs  = 30

tmoFiles = ['/Users/biggus/Documents/James/Data/ReClustering/PrelimData_Grant_Feb09/RandSplitFastas/MemeResults/Clus2_247gene_0.8_Apr16_14-46-36.meme.txt.tmo',
            '/Users/biggus/Documents/James/Data/ReClustering/PrelimData_Grant_Feb09/RandSplitFastas/MemeResults/Clus2_247gene_0.8_Apr16_14-46-33.meme.txt.tmo',]

fastaPath   = '/Users/biggus/Documents/James/Data/2KB/2kb_Sequence/2kb_Anopheles/2KBupTSS_goodAffyAGAPsFastasOUT.masked.nr.fas'
seqs        = Fasta.load(fastaPath)

targetGenes = '/Users/biggus/Documents/James/Data/ReClustering/kmedsPear33Clus50x_2/Clus2_247genes.genes.txt'
targetGenes = map(lambda l: l.strip(), open(targetGenes, 'rU'))
targetGenes = targetGenes[:genes]
for i in range(len(targetGenes)):
    targetGenes[i] = seqs[targetGenes[i]]


motifs  = []
tMotifs = []
mMotifs = []
for t in tmoFiles:
    Ms = loadTMOs(t)
    motifs.extend(Ms)
for i in range(len(motifs)):
    tMotifs.append(Motif(motifs[i].bogus_kmers()))
    mMotifs.append(motility.make_pwm(motifs[i].bogus_kmers()))



tTimeIt =\
"""for m in tMotifs:
    for s in targetGenes:
        m.scan(s,factor=0.75)
"""
mTimeIt =\
"""for m in mMotifs:
    for s in targetGenes:
        m.find(s,m.max_score()*0.85)
"""


tTimer = timeit.Timer(tTimeIt,"from __main__ import tMotifs,targetGenes")
mTimer = timeit.Timer(mTimeIt,"from __main__ import mMotifs,targetGenes")

tTimes = tTimer.repeat(repeat=runs,number=1)
mTimes = mTimer.repeat(repeat=runs,number=1)

tOut = open(tOut, 'w')
mOut = open(mOut, 'w')

for t in tTimes:
    tOut.write('%.6f\n' % (t))
for t in mTimes:
    mOut.write('%.6f\n' % (t))
    
print "IM DONE!!!  "*20
        
        