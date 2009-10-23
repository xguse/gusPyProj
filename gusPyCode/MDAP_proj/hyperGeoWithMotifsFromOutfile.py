# get hypergeo pVal from motifs from MD output files

from gusPyCode.MDAP_proj.MDAP_defs import loadMotifsFromOutFile
from TAMO.MotifMetrics import ProbeSet
from gusPyCode.defs.mosqData import promoterSeqPaths

outFile   = '/Users/biggus/Documents/James/Data/ReClustering/PrelimData_Grant_Feb09/Clus2_kmerSearch.6-8mers.FDR_lessThan0.01.analysis.txt'

fileNames = ['/Users/biggus/Documents/James/Data/ReClustering/PrelimData_Grant_Feb09/Clus2_kmerSearch.6-8mers.FDR_lessThan0.01.txt',]

coRegSeqs    = '/Users/biggus/Documents/James/Data/ReClustering/kmedsPear33Clus50x_2/Clus2_247genes.genes.txt'
allPromoters = promoterSeqPaths.Aa_2000bpUp_hardMasked


MDfileType = 'list'
dfltFactor = 0.75

coRegSeqs = map(lambda l: l.strip(), open(coRegSeqs, 'rU').readlines())

motifs = []

for f in fileNames:
    motifs.append([f.split('/')[-1],loadMotifsFromOutFile(f, MDfileType)])

probSet = ProbeSet(allPromoters)

# get and print motif and pVal(s)
out = '#outFile\tMotif\tHyperGeoPval (%s)\tfrac (%s)\tBestHG (pval)\tfrac (bestHG_ScoreThresh)\tBestHG (scoreThresh)\tbinoPval (%s)\tbinoPval (bestHG_scoreThresh)' \
      % (dfltFactor,dfltFactor,dfltFactor)
print out
out = out+'\n'

for m in motifs:
    for i in m[1]:
        bestE = probSet.best_p_value(m[1],coRegSeqs)
        temp ='%s\t%s\t%.3e\t%.3f\t%.3e\t%.3f\t%.3f\t%.3e\t%.3e' % (m[0],
                                                      m[1], 
                                                      probSet.Enrichment(m[1],coRegSeqs,factor=dfltFactor), 
                                                      probSet.frac(m[1],coRegSeqs,factor=dfltFactor),
                                                      bestE[0],
                                                      probSet.frac(m[1],coRegSeqs,factor=bestE[1]),
                                                      bestE[1],
                                                      probSet.binomial(m[1],coRegSeqs,factor=dfltFactor),
                                                      probSet.binomial(m[1],coRegSeqs,factor=bestE[1]))
        #bestE = probSet.best_p_value(i,coRegSeqs)
        #temp ='%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (m[0],
                                                      #i, 
                                                      #probSet.Enrichment(i,coRegSeqs,factor=dfltFactor), 
                                                      #probSet.frac(i,coRegSeqs,factor=dfltFactor),
                                                      #bestE[0],
                                                      #probSet.frac(i,coRegSeqs,factor=bestE[1]),
                                                      #bestE[1],
                                                      #probSet.binomial(i,coRegSeqs,factor=dfltFactor),
                                                      #probSet.binomial(i,coRegSeqs,factor=bestE[1]))
        print temp
        out = out+temp+'\n'
        

outFile = open(outFile,'w')
outFile.write(out)
outFile.close()

print 'done'
