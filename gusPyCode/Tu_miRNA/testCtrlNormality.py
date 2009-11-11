from random import choice
from scipy.stats import normaltest
from numpy import log10, mean, median
import cPickle
from gusPyCode.defs import miRNA_targeting as miRT


print '\n\n'

outFile = '/Users/biggus/Documents/James/Data/Tu_miRNA/ResultsMiRNA_Targeting/2009_08_31run.100ctrls.ctrlNormality.AP.txt'
outFile = open(outFile, 'w')

pklPath = '/Users/biggus/Documents/James/Data/Tu_miRNA/ResultsMiRNA_Targeting/2009_08_31.seedMatches.100ctrls.pkl'
data    = cPickle.load(open(pklPath,'rU'))

# calc the various centers and spreads
#items = ['aga-miR-1890',
         #'aga-miR-9c',
         #'aga-miR-13b',
         #'aga-miR-133',
         #'aga-miR-11',
         #'aga-miR-988',
         #'aga-miR-278',
         #'aga-miR-965',]

#miR_matches = {}
#for i in items:
    #miR_matches[i] = data['miR_matches'][i]

miR_matches = data['miR_matches']
    
def writeResults(miRNA,outFile,log=0):
    """
    write out each miRNA seedType's normaltest results for each orthoType.
    """
    print miRNA.name
    for seedType in sorted(miRNA.matchVersions.keys()):
        for orthoType in range(4):
            data = []
            if log: data.extend([log10(x[orthoType]) for x in miRNA.ctrlCounts[seedType]])
            else: data.extend([x[orthoType] for x in miRNA.ctrlCounts[seedType]])
            results = normaltest(data)
            outFile.write('%s\t%s\t%s\t%s\t%.2F\t%.2F\t%.5g\t%.5g\n' % \
                          (miRNA.name,
                           seedType,
                           orthoType,
                           len(data),
                           mean(data),
                           median(data),
                           results[0],
                           results[1]))
        outFile.flush()
            

outFile.write('miRNA\tseedType\tOrthoType\tN\tmean\tmedian\tW-stat\tpVal4W\n')
for miRNA in sorted(miR_matches.keys()):
    writeResults(miR_matches[miRNA],outFile,log=0)