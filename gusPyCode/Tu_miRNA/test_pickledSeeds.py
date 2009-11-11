import cPickle
from matplotlib import pylab as plb
from gusPyCode.defs import miRNA_targeting as miRT

print '\n\n'

outFile = '/Users/biggus/Documents/James/Data/Tu_miRNA/ResultsMiRNA_Targeting/2009_09_29run.ps100ctrls.compStats.txt' #'/Users/biggus/Documents/James/Data/Tu_miRNA/ResultsMiRNA_Targeting/2009_08_31run.100ctrls.stats.txt'
outFile = open(outFile, 'w')

pklPath = '/Users/biggus/Documents/James/Data/Tu_miRNA/ResultsMiRNA_Targeting/2009_09_29.seedMatches.100psCtrls.pkl' #'/Users/biggus/Documents/James/Data/Tu_miRNA/ResultsMiRNA_Targeting/2009_08_31.seedMatches.100ctrls.pkl'
data    = cPickle.load(open(pklPath,'rU'))

# calc the various centers and spreads

seeds = data['miR_matches']
print "Getting centers and spreads..."
for seed in seeds:
    seeds[seed].calcCtrlMeanStDv()
    seeds[seed].calcCtrlMedianStDv()
    seeds[seed].calcCtrlMedianMAD()
    
# Write out file documenting each seed's seedType stats
print "Writting stats for miRNA:"
c = 0
seedTypes = sorted(miRT._seedModels.keys())
outFile.write('#miRNA\tseedType\torthoType\trealData\tmean\tmedian\tmeanStDv\tmedianStDv\tmedianMAD\tmeanStDvZ-score\tmedianStDvZ-score\tmedianMADZ-score\n')
for seed in sorted(seeds.keys()):
    c+=1
    print c
    for seedType in seedTypes:
        for orthoNum in range(4):
            outFile.write('%s\t%s\t%s\t%s\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n' \
                          % (seed,
                             seedType,
                             orthoNum,
                             seeds[seed].matchCounts[seedType][orthoNum],
                             seeds[seed].ctrlMeanStd[seedType][orthoNum][0],
                             seeds[seed].ctrlMedianStd[seedType][orthoNum][0],
                             seeds[seed].ctrlMeanStd[seedType][orthoNum][1],
                             seeds[seed].ctrlMedianStd[seedType][orthoNum][1],
                             seeds[seed].ctrlMedianMAD[seedType][orthoNum][1],
                             seeds[seed].scoreSeedType(seedType, orthoNum,'meanStDv'),
                             seeds[seed].scoreSeedType(seedType, orthoNum,'medianStDv'),
                             seeds[seed].scoreSeedType(seedType, orthoNum,'medianMAD')))
            
outFile.flush()
outFile.close()
print 'Im Done.'