import cPickle, time
import numpy
from TAMO.MotifTools import Motif
from TAMO import MotifTools
from TAMO.Clustering import MotifCompare
from gusPyCode.MDAP_proj.MDAP_defs import getMinDiffOri,getKmersWithOneMisMtch,alignSimilarMotifs,trimPaddedMotif
## NOTE: defs moved to MDAP_defs
#def getKmersWithOneMisMtch(motif1, motifListWithMetrics):
    #"""
    #Takes a TAMO motif and a list of lists: [TAMOmotif, weightMetric].
    #Returns listOfLists in same form but containing only those motifs of
    #length = len(motif1), and if a revComp in motifListWithMetrics matches
    #better, _IT_ is returned instead of the original motif.
    #"""
    #resultList = []
    
    #for mWithMetric in motifListWithMetrics:
        ## Determine what distanceResult == one misMatch.
        ## Use length of shortest motif for misMatch Calc
        #maxAlnLen   = min(len(motif1), len(mWithMetric[0]))
        #kMer        = Motif('%s' %('A'*maxAlnLen))
        #kMer1mis    = Motif('%s%s' %('A'*(maxAlnLen-1),'T'))
        #oneMisMatch = MotifCompare.minshortestoverhangdiff(kMer,kMer1mis)
        #bestOri = getMinDiffOri(motif1,mWithMetric[0])
        #if bestOri[1] <= oneMisMatch:
            #resultList.append([bestOri[0],mWithMetric[1]]) # keep [bestOriTAMOmotif, weightMetric]
            
    #return resultList

#def alignSimilarMotifs(motifsWithOneMisMtch):
    #"""
    #Takes list output from getKmersWithOneMisMtch(). List should be sorted by weight.
    #Pads top motif excessivly, then pads all other motifs to align with top motif.
    #Resturns list of padded, aligned motifs.
    #"""
    ## Pad top motif with 5 on the left side, 5 on the right.
    #motifsWithOneMisMtch[0][0] = motifsWithOneMisMtch[0][0][-5,len(motifsWithOneMisMtch[0][0])+5]
    
    #for i in range(1,len(motifsWithOneMisMtch)):
        #offSet = MotifCompare.minshortestoverhangdiff(motifsWithOneMisMtch[0][0],motifsWithOneMisMtch[i][0], want_offset=1)
        ## Pad to align to top motif
        #motifsWithOneMisMtch[i][0] = motifsWithOneMisMtch[i][0][-offSet[0],len(motifsWithOneMisMtch[i][0])]
        #motifsWithOneMisMtch[i][0] = motifsWithOneMisMtch[i][0][0,len(motifsWithOneMisMtch[0][0])]
    #return motifsWithOneMisMtch

#def trimPaddedMotif(paddedMotif):  # doesnt seem needed since motifTools.sum() seems to do this automatically
    #"""
    #Takes a padded motif and removes blank positions from left and right, until it
    #sees a non-blank.  Returns trimmed motif obj.
    #"""
    
    ## Count Blanks
    #lBlanks     = 0
    #rBlanks     = 0
    #posInfoBits = paddedMotif.bits[:] # real copy bc i will be playing with orders and dont want to cause trouble
    #for pos in posInfoBits:
        #if pos == 0:
            #lBlanks+=1
        #else:
            #break
    #posInfoBits.reverse()   
    #for pos in posInfoBits:
        #if pos == 0:
            #rBlanks+=1
        #else:
            #break
        
    ## Remove Blanks
    #motif = paddedMotif[lBlanks:-rBlanks+1]
    
    #return motif
    
        

# # # # # # # 
# code for testing defs above

#mo = Motif(Motif('....WGATAAR').bogus_kmers())
#print mo.print_textlogo()

t1 = time.time()
testMotifs = '/Users/biggus/Documents/James/Collaborations/Campbell/data/Results_HyperGeoScreen/unMasked/CCupAt4Days.genes.6-9mers.UnMasked.txt'
testMotifs = map(lambda l: l.strip().split('\t'), open(testMotifs, 'rU').readlines())

if testMotifs[0][0].startswith('#'): testMotifs.pop(0) # remove header if present

# TAMOify kmers and logify pVals
for i in range(len(testMotifs)):
    testMotifs[i] = (Motif(testMotifs[i][0]),numpy.log10(float(testMotifs[i][1])))
    
# Sort on log'd pVals
testMotifs.sort(key=lambda x: x[1])

comboMotifs = []

for i in range(0,int(len(testMotifs)*0.2)):
    simMotifs = getKmersWithOneMisMtch(testMotifs[i][0],testMotifs) 
    simMotifs = alignSimilarMotifs([x[0] for x in simMotifs])
    #for m in simMotifs:
        #print m[0].oneletter
    comboMotifs.append(MotifTools.sum([x[0] for x in simMotifs],[-x[1] for x in simMotifs])) # -x[1] to convert neg logs to pos weights
    print len(comboMotifs)

t2 = time.time()    

oFile = '/Users/biggus/Desktop/CCupAt4Days.genes.6-9mers.UnMasked.gGEMS.tmo'
pFile = '/Users/biggus/Desktop/CCupAt4Days.genes.6-9mers.UnMasked.gGEMS.pkl'
MotifTools.save_motifs(comboMotifs,oFile,kmer_count=60)

pFile = open(pFile, 'w')
cPickle.dump(comboMotifs,pFile)
t3 = time.time()    
print 'Calculations took %.3f min.\nWriting/Pickling took %.3f min.' % ((float(t2)-t1)/60, (float(t3)-t2)/60) 
    


    
None


