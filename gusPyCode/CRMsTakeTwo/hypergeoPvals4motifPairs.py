from time import time
import random
import cPickle
from TAMO.MotifTools import load as loadTMOs
from crmClasses import *

mapPickle = '/Users/biggus/Desktop/testMap.pkl'
pValOut   = '/Users/biggus/Desktop/testPvals.txt'

fastaPath   = '/Users/biggus/Documents/James/Data/2KB/2kb_Sequence/2kb_Anopheles/2KBupTSS_goodAffyAGAPsFastasOUT.masked.nr.fas'
targetGenes = '/Users/biggus/Documents/James/Data/ReClustering/kmedsPear33Clus50x_2/Clus2_247genes.genes.txt'

tPkls    = None # ['/Users/biggus/Documents/James/Data/JasparMotifs/JASPAR_CORE_2008/jsprCore2008.tamo.pkl',]

tmoFiles = ['/Users/biggus/Documents/James/Data/ReClustering/PrelimData_Grant_Feb09/RandSplitFastas/AceResults/Clus2_247gene_0.8_Apr16_14-46-36.ace.2.txt.tmo',
            '/Users/biggus/Documents/James/Data/ReClustering/PrelimData_Grant_Feb09/RandSplitFastas/MemeResults/Clus2_247gene_0.8_Apr16_14-46-36.meme.txt.tmo',
            '/Users/biggus/Documents/James/Data/ReClustering/PrelimData_Grant_Feb09/RandSplitFastas/MemeResults/Clus2_247gene_0.8_Apr16_14-46-33.meme.txt.tmo',
            '/Users/biggus/Documents/James/Data/ReClustering/PrelimData_Grant_Feb09/RandSplitFastas/AceResults/Clus2_247gene_0.8_Apr16_14-46-36.ace.3.txt.tmo',
            '/Users/biggus/Documents/James/Data/ReClustering/PrelimData_Grant_Feb09/RandSplitFastas/AceResults/Clus2_247gene_0.8_Apr16_14-46-36.ace.1.txt.tmo',
            '/Users/biggus/Documents/James/Data/ReClustering/PrelimData_Grant_Feb09/RandSplitFastas/AceResults/Clus2_247gene_0.8_Apr16_14-46-33.ace.3.txt.tmo',
            '/Users/biggus/Documents/James/Data/ReClustering/PrelimData_Grant_Feb09/RandSplitFastas/AceResults/Clus2_247gene_0.8_Apr16_14-46-33.ace.2.txt.tmo',
            '/Users/biggus/Documents/James/Data/ReClustering/PrelimData_Grant_Feb09/RandSplitFastas/AceResults/Clus2_247gene_0.8_Apr16_14-46-33.ace.1.txt.tmo',
            '/Users/biggus/Documents/James/Data/ReClustering/PrelimData_Grant_Feb09/Clus2_247genes.6-8mers.gGEMS.tmo']


targetGenes = map(lambda l: l.strip(), open(targetGenes, 'rU'))


lenFiles = 0
if tPkls: lenFiles += len(tPkls)
if tmoFiles: lenFiles += len(tmoFiles)



motifs = []
if tPkls != None: 
    for p in tPkls:
        Ms = cPickle.load(open(p, 'r'))
        for i in range(len(Ms)):
            Ms[i].sourceFile = p
        motifs.extend(Ms)
        print '%s motifs from %s' % (len(Ms), p.split('/')[-1])
if tmoFiles != None:
    for t in tmoFiles:
        Ms = loadTMOs(t)
        for i in range(len(Ms)):
            Ms[i].sourceFile = t
        motifs.extend(Ms)
        print '%s motifs from %s' % (len(Ms), t.split('/')[-1])

# to debug run time/motif         
motifDict = {}
mCount = 0
for m in motifs:
    mCount+=1
    motifDict['%03d_%s' % (mCount, m.oneletter)] = m
    
#rMos = []
#for i in range(16):
    #rMos.append(random.choice(motifDict.keys()))
#motifs = rMos
#for m in motifs:
    #print m

motifs = ['466_.aACCCA.',
          '395_.TATCAr',
          '120_Gs.G...G..sCG...G..sG',
          '508_cCAMCCg',
          '071_TG.TG.TG.TG.TG',
          '153_GGGrrGrrGG',
          '198_AA.sAAA..AA....AA',
          '086_r.r.GGGGAGGG',
          '374_.tcwtat',
          '498_.cCGATc',
          '146_g.r.GGGGAGGg',
          '069_mAAACAAmAm',
          '482_CGATC.',
          '258_CSCGc..GCG.kcG',
          '397_..CGAKm',
          '339_Gs..Gs.GC.kGC.G',]   

tmpDict = {}
for m in motifs:
    tmpDict[m] = motifDict[m]
motifDict = tmpDict
# end debug code


print '%s motifs loaded from %s files.' % (len(motifDict), lenFiles)
t1 = time()
maps = MapLib(fastaPath,motifDict, thresh=0.7)
t2 = time()
cPickle.dump(maps,open(mapPickle, 'w'))
print 'Map creation took %.2f min.' % ((t2-t1)/60.0)
t3 = time()
tGeneMotPairs = maps.getMotifPairs4GeneList(gList)
t4 = time()
print 'Generating target motif pairs took %.2f min.' % ((t4-t3)/60.0)
pairPvals = {}
t5 = time()
for pair in tGeneMotPairs:
    pairPvals[pair] = maps.calcPairsCumHG_pVal(gList,pair)
t6 = time()
print 'Calculating p-values took %.2f min.' % ((t4-t3)/60.0)


    
    