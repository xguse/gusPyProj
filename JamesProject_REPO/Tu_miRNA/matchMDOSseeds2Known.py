from JamesDefs import loadXXmiRNAs
from JamesDefs import revComp


kMersFile     = '/Users/biggus/Documents/James/Data/Tu_miRNA/MDOSoutPut/CqAg_7mers_orthosForMDOS_500afterCoding.mdos.motifzSgte3.txt'
knownRNAsPATH = '/Users/biggus/Documents/James/Data/Tu_miRNA/miRNAs/Ag_miRNAa_unknSource/aga.fa'
oFile         = '/Users/biggus/Documents/James/Data/Tu_miRNA/MDOSoutPut/MDOS_kmer_IDs/CqAg_7mers_vs_AgMiRNAs.ctrl.txt'

knownRNAs= loadXXmiRNAs(knownRNAsPATH)
kMers = map(lambda line: line.strip().split('\t'), open(kMersFile,'rU').readlines())
mdosFieldNames = kMers.pop(0)

oList = ['# kMersFile: %s\n# knownRNAsPATH: %s\n\n%s\n' % (kMersFile,knownRNAsPATH,mdosFieldNames)]



for k in kMers:
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
    else:
        oList.append('NO MATCH\n')
    

            
oFile = open(oFile, 'w')
oFile.writelines(oList)

print 'done'