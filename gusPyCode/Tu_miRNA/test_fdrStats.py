import cPickle
from gusPyCode.defs import miRNA_targeting as miRT
from gusPyCode.defs.mathDefs import mean

print '\n\n'

#outFile = '/Users/biggus/Documents/James/Data/Tu_miRNA/ResultsMiRNA_Targeting/testPickle.5randMiRs.100prmCtrls.stats.events.all.txt'
#outFile = open(outFile, 'w')

pklPath = '/Users/biggus/Documents/James/Data/Tu_miRNA/ResultsMiRNA_Targeting/2009_08_31.seedMatches.100ctrls.pkl'
data    = cPickle.load(open(pklPath,'rU'))

miR_matches = data['miR_matches']
for miR in miR_matches:
    print miR
    for seedType in miRT._seedModels:
        for i in range(1,4):
            sts = miR_matches[miR].getFDRStats(seedType,i)
            realMatches = miR_matches[miR].matchCounts[seedType][i]
            if sts[0] != None:
                if (sts[0]+(3*sts[1])) < 0.15:
                    print '\t%s[%s]:%.3G - %s hits' % (seedType,i,sts[0],realMatches)
            
    