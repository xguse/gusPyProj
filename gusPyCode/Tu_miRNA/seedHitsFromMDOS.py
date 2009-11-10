from TAMO.MotifMetrics import Fasta
from gusPyCode.defs.JamesDefs import revComp

seedsFile = '/Users/biggus/Documents/James/Data/Tu_miRNA/miRNAs/testSeeds.fas'
kMersFile = '/Users/biggus/Documents/James/Data/Tu_miRNA/MDOSoutPut/CqAg_7mers_orthosForMDOS_500afterCoding.mdos.motifzSgte3.txt'

seeds = Fasta.file2dict(seedsFile)
kMers = map(lambda line: line.strip().split('\t'), open(kMersFile,'rU').readlines())

print 'seedsFile: %s\nkMersFile: %s\n' % (seedsFile,kMersFile)

for l in kMers:
    for seed in seeds.keys():
        if revComp(seeds[seed].upper()) == l[0]:
            print '>%s\n%s\n' % (seed,'\t'.join(l))

            
print 'done'