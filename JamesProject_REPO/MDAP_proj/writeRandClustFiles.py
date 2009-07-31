from MDAP_defs import genRandClusters
from TAMO.seq import Fasta
from numpy import average
import pprint
import glob



clusterFile = '/Users/biggus/Documents/James/Data/ReClustering/kmedsPear33Clus50x_2/Clus2_247genes.genes.txt'
totalSeqs = '/Users/biggus/Documents/James/Data/2KB/2kb_Sequence/2kb_Anopheles/2KBupTSS_goodAffyAGAPsFastasOUT.masked.nr.fas'
totalSeqs = Fasta.file2dict(totalSeqs)

geneNames = map(lambda l: l.strip(),open(clusterFile, 'rU').readlines())
randClusterLists = genRandClusters(geneNames,totalSeqs,N=3, keepLen=1)

for i in range(len(randClusterLists)):
    oFile = open(clusterFile.replace('genes.txt','rGenes%s.txt' % (i)), 'w')
    for name in randClusterLists[i]:
        oFile.write(name+'\n')
    oFile.close()
    del(oFile)
    
print 'Done.'


            
print "Original list:"
origLens = []
for g in geneNames:
    origLens.append(len(totalSeqs[g]))
origLens.sort()
print origLens
print 'Avg: %.3f' % (average(origLens))
    
print "Randomized lists:"
for r in randClusterLists:
    randLens = []
    for g in r:
        randLens.append(len(totalSeqs[g]))
    randLens.sort()
    print randLens
    print 'Avg: %.3f' % (average(randLens))
    