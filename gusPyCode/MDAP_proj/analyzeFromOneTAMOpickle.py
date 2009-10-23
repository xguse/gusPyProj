# get hypergeo pVal from motifs from TAMO motifs pickles

from MDAP_defs import loadMotifsFromOutFile
from TAMO.MotifMetrics import ProbeSet
import cPickle
import pprint

outFile   = '/Users/biggus/Documents/James/Data/ReClustering/PrelimData_Grant_Feb09/Clus2_247.kmerSearch.gGEM.analysis.txt'

pklFilePath = '/Users/biggus/Documents/James/Data/ReClustering/PrelimData_Grant_Feb09/Clus2_247.kmerSearch.pkl'
pklFile = open (pklFilePath, 'r')

coRegSeqs = '/Users/biggus/Documents/James/Data/ReClustering/kmedsPear33Clus50x_2/Clus2_247genes.genes.txt'

dfltFactor = 0.75

coRegSeqs = map(lambda l: l.strip(), open(coRegSeqs, 'rU').readlines())

motifs = cPickle.load(pklFile) 
for i in range(len(motifs)):
    motifs[i] = [pklFilePath.split('/')[-1],motifs[i]]
    



probSet = ProbeSet('/Users/biggus/Documents/James/Data/2KB/2kb_Sequence/2kb_Anopheles/2KBupTSS_goodAffyAGAPsFastasOUT.masked.nr.fas')

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
