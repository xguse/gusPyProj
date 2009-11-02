print 'loading modules...'
from gusPyCode.gClasses.gSeqClasses import DNAseq
import motility
from matplotlib import pyplot
from pprint import pprint



print 'reading in fasta file...'
alignedSeqFile  = map(lambda line: line.strip('\n') , open('/Users/biggus/Documents/James/Data/AedesCorePromoterWork/output/TATA_Inr/Aa_Ensbl49_AaegL1.1.plus100minus100._TataInr_.newInr.geneStrand.fa','rU').readlines())
alignedSeqs     = zip(alignedSeqFile[:-1:2], alignedSeqFile[1::2])
motif1 = motility.IUPAC('TATAWAAR')
motif2 = motility.IUPAC('TCAKTY')
motif1.name,motif2.name = 'TATA','Inr'
xAxisText = 'Nucleotide position'
yAxisText = 'Motif occurences'

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
yData1 = []
yData2 = []
for i in range(0,max(seqLengths)):
    yData1.append(0)
    yData2.append(0)

for name in alignedSeqsDict:
    hits = motif1.find(alignedSeqsDict[name].toString())
    for hit in hits:
        if hit[2] == 1:
            yData1[int(hit[0])] += 1
            
for name in alignedSeqsDict:
    hits = motif2.find(alignedSeqsDict[name].toString())
    for hit in hits:
        if hit[2] == 1:
            yData2[int(hit[0])] += 1

pyplot.xlabel(xAxisText)
pyplot.ylabel(yAxisText)


pyplot.plot(xData,yData1,'b',label=motif1.name)
pyplot.plot(xData,yData2,'r',label=motif2.name)
pyplot.legend(loc='upper right')


pyplot.show()

