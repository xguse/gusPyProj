# get hypergeo pVal from motifs from TAMO motifs pickles
from gusPyCode.defs import nucFreqRepo
from gusPyCode.MDAP_proj.MDAP_defs import loadMotifsFromOutFile
from TAMO.MotifMetrics import ProbeSet
import cPickle
import pprint

outFile   = '/Users/biggus/Documents/James/Collaborations/Campbell/data/Results_HyperGeoScreen/masked/Results_gGEMS/CCupAt4Days.gte2x.5-16mers.gGEMS.analysis70.txt'

pklFilePath = '/Users/biggus/Documents/James/Collaborations/Campbell/data/Results_HyperGeoScreen/masked/Results_gGEMS/CCupAt4Days.gte2x.5-16mers.gGEMS.pkl'
pklFile = open (pklFilePath, 'r')

coRegSeqs = '/Users/biggus/Documents/James/Collaborations/Campbell/data/CCupAt4Days.gte2x.genes.txt'

dfltFactor = 0.70
speciesBK  = nucFreqRepo.AaMasked_NucFreqs_2kUpPromo

coRegSeqs = map(lambda l: l.strip(), open(coRegSeqs, 'rU').readlines())

motifs = cPickle.load(pklFile) 
for i in range(len(motifs)):
    motifs[i] = [pklFilePath.split('/')[-1],motifs[i]]
    
# Adjust each motif to the species being looked at
for i in motifs:
    motifs[i][1].new_bg(speciesBK)



probSet = ProbeSet('/Users/biggus/Documents/James/Data/2KB/2kb_Sequence/2kb_Aedes/aedes2KBupStreamTSS.hardMasked.geneStrand.fas')

# get and print motif and pVal(s)
out = '#outFile\tMotif\tHyperGeoPval (%s)\tfrac (%s)\tBestHG (pval)\tfrac (bestHG_ScoreThresh)\tBestHG (scoreThresh)\tbinoPval (%s)\tbinoPval (bestHG_scoreThresh)' \
      % (dfltFactor,dfltFactor,dfltFactor)
print out
out = out+'\n'


for m in motifs:
    bestE = probSet.best_p_value(m[1],coRegSeqs)
    temp ='%s\t%s\t%.3e\t%.3f\t%.3e\t%.3f\t%.3f\t%.3e\t%.3e' % (m[0],
                                                  m[1].oneletter, 
                                                  probSet.Enrichment(m[1],coRegSeqs,factor=dfltFactor), 
                                                  probSet.frac(m[1],coRegSeqs,factor=dfltFactor),
                                                  bestE[0],
                                                  probSet.frac(m[1],coRegSeqs,factor=bestE[1]),
                                                  bestE[1],
                                                  probSet.binomial(m[1],coRegSeqs,factor=dfltFactor),
                                                  probSet.binomial(m[1],coRegSeqs,factor=bestE[1]))
    print temp
    out = out+temp+'\n'
    

outFile = open(outFile,'w')
outFile.write(out)
outFile.close()

print 'done'
