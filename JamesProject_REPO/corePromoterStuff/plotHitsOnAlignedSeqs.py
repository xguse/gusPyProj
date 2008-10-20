print 'loading modules...'
from gSeqClasses import DNAseq
import motility
from matplotlib import pylab
from pprint import pprint

print 'reading in fasta file...'
alignedSeqFile  = map(lambda line: line.strip('\n') , open('/Users/biggus/Documents/James/Data/AedesCorePromoterWork/output/Aa_Ensbl49_AaegL1.1.plus50minus100._InrDPE_seqSlices_.fa','rU').readlines())
alignedSeqs     = zip(alignedSeqFile[:-1:2], alignedSeqFile[1::2])
motif = motility.IUPAC('RGWYV')

alignedSeqsDict = {}
# populate alignedSeqsDict
print 'populating alignedSeqsDict...'
for item in alignedSeqs:
    name = item[0].replace('>','')
    seq  = item[1]
    alignedSeqsDict[name] = DNAseq(seq,name)
print ' ...alignedSeqsDict has %s entries.' % (len(alignedSeqsDict))

seqLengths = []
for s in alignedSeqsDict:
    seqLengths.append(len(alignedSeqsDict[s]))
    
xData = range(0,max(seqLengths))
yData = []
for i in range(0,max(seqLengths)):
    yData.append(0)

for name in alignedSeqsDict:
    hits = motif.find(alignedSeqsDict[name].toString())
    for hit in hits:
        if hit[2] == 1:
            yData[int(hit[0])] += 1

pylab.xlabel('position relative to end of Inr motif')
pylab.ylabel('DPE motifs')

pylab.plot(xData,yData,'b')

pylab.show()

