print "loading modules..."
from TAMO import MotifTools 
from TAMO.seq import Fasta 
from TAMO import MotifMetrics
from TAMO.MD.AlignAce import AlignAce 
from TAMO.MD.MDscan import MDscan 
from TAMO.MD.Meme import Meme 
from TAMO.Clustering import Kmedoids,MotifCompare,UPGMA
#from TAMO.DataSources import GO
from time import time
from pprint import pprint

print "loading vars..."
motifList = map(lambda line: MotifTools.Motif_from_text(line.strip()),open(\
    '/Users/biggus/Documents/James/Writings_Talks/Grants/09_Feb/PrelimData_Grant_Feb09/Clus2_kmerSearch-0.01.8mers.motifs.txt','r').readlines())

dMat      = ''

print "constructing distanceMatrix..."
dM_t1 = time()
distanceMatrix={}
for i in range(len(motifList)):
    print 'motif %s of %s' % (i+1,len(motifList))
    distanceMatrix[i]={}
    for j in range(len(motifList)):
	# check fwd and revCmp alignments and take the lowest
	fwd_diff    = MotifCompare.minshortestoverhangdiff(motifList[i],motifList[j])
	revCmp_diff = MotifCompare.minshortestoverhangdiff(motifList[i].revcomp(),motifList[j])
	distanceMatrix[i][j] = min([fwd_diff,revCmp_diff])
dM_t2 = time()
pprint(distanceMatrix)
print 'distanceMatrix took %.4f sec.' % (dM_t2-dM_t1) 
		
print "discovering clusters..."
# --Using Kmedoids --
clusterOut = Kmedoids.bestaveKMedoids_cluster(distanceMatrix,kmax=30)

for c in clusterOut[1]:
    print 'cluster_%s:' % (c)
    for m in clusterOut[1][c]:
	print motifList[m].oneletter
    print '\n'

## --Using UPGMA --
#print "Clustering %d motifs"%len(motifList)
#clsr_t1 = time()
#DMAX  = 0.22

#tree = UPGMA.UPGMA(motifList)
#print 'printing tree:'
#UPGMA.print_tree(tree)


#motifs = UPGMA.slice_tree(tree,DMAX)
#for m in motifs:
    #print m.oneletter,m.revcomp().oneletter,len(UPGMA.flatten_members(m))
#clsr_t2 = time()
#print 'Clustering took %.3f sec.' % (clsr_t2-clsr_t1)



x=1