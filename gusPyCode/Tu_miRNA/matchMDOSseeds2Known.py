from gusPyCode.defs.JamesDefs import loadXXmiRNAs
from gusPyCode.defs.JamesDefs import revComp


kMersFile     = '/Users/biggus/Documents/James/Data/Tu_miRNA/MDOSoutPut/AaAgAndAgCq_7mers_orthosForMDOS_500afterCoding.mdos.allMotifsNoAmbig.txt'
knownRNAsPATH = '/Users/biggus/Documents/James/Data/Tu_miRNA/miRNAs/Ag_miRNAa_unknSource/aga.fa'
oFile         = '/Users/biggus/Documents/James/Data/Tu_miRNA/MDOSoutPut/MDOS_kmer_IDs/AaAgAndAgCq_7mers_vs_AgMiRNAs.allMotifsNoAmbig.txt'
histDataFile  = '/Users/biggus/Documents/James/Data/Tu_miRNA/MDOSoutPut/MDOS_kmer_IDs/AaAgAndAgCq_7mers_vs_AgMiRNAs.allMotifsNoAmbig.histData.txt'

knownRNAs= loadXXmiRNAs(knownRNAsPATH)
kMers = map(lambda line: line.strip().split('\t'), open(kMersFile,'rU').readlines())
mdosFieldNames = kMers.pop(0)

oList = ['# kMersFile: %s\n# knownRNAsPATH: %s\n\n%s\n' % (kMersFile,knownRNAsPATH,mdosFieldNames)]



KmerZscores = []
hitZscores  = []
missZscores = []

for k in kMers:
    KmerZscores.append(k[7])
    oList.append('--\n%s\n' % ('\t'.join(k)))
    kHits = []
    
    for r in knownRNAs:
        pos2 = r[2][-8:-1] # represents hit to seed starting at pos2
        pos1 = r[2][-7:]   # hit starting at pos1
        
        if k[0] == pos2:   # to search ctrl use ''' if revComp(k[0]) == pos2: '''
            pos2 = 'yes'
        else:
            pos2 = 'no'
            
        if k[0] == pos1:   # to search ctrl use ''' if revComp(k[0]) == pos1: '''
            pos1 = 'yes'
        else:
            pos1 = 'no'
        
        if pos1 == 'yes' or pos2 == 'yes':
            kHits.append('>%s\n' % (r[0]))
            kHits.append('%s\tmiRNA\n' % (r[1]))
            kHits.append('%s\trevComp\tpos2(%s)\tpos1(%s)\n\n' % (r[2],pos2,pos1))
            
    if kHits:
        oList.extend(kHits)
        hitZscores.append(k[7])
    else:
        oList.append('NO MATCH\n')
        missZscores.append(k[7])
        
    

assert len(hitZscores)+len(missZscores) == len(KmerZscores), \
       '''len(hitZscores)+len(missZscores) != len(KmerZscores)
       hitZscores  = %s
       missZscores = %s
       KmerZscores = %s''' \
       % (len(hitZscores),len(missZscores),len(KmerZscores))
        
oFile = open(oFile, 'w')
oFile.writelines(oList)

histDataFile = open(histDataFile,'w')
histDataFile.write('\t'.join(KmerZscores)+'\n')
histDataFile.write('\t'.join(hitZscores)+'\n')

print 'Total kmers: %s.\nKmers matching at least one DB entry: %s.' % (len(kMers),len(hitZscores))

print 'done'