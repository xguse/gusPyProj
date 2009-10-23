print 'loding modules...'
import re
from gSeqClasses import DNAseq
from gusPyCode.defs.JamesDefs import overlapRegEx
from gusPyCode.defs.bioDefs import iupac2regex

print 'loading fastaFile...'
seqFile     = map(lambda line: line.strip(), open('/Users/biggus/Documents/James/Data/AedesCorePromoterWork/output/Aa_Ensbl49_AaegL1.1.plus50minus100._InrDPE_Hits_.fas','rU').readlines())
seqs        = zip(seqFile[:-1:2], seqFile[1::2])
motifCombo  = re.compile(iupac2regex('TCAKTYN{23}RGWYV'),re.I)

seqsDict = {}
# populate seqsDict
print 'populating seqsDict...'
for item in seqs:
    name   = item[0].split(' ')[0].split(':')[0].replace('>','')
    header = item[0]
    seq    = item[1]
    seqsDict[name] = DNAseq(seq,name)
print 'loaded %s entries into seqsDict.' % (len(seqsDict))

# finds non-overlapping patters
##hitCount = {}
##for name in seqsDict:
    ##hits = motifCombo.finditer(seqsDict[name].toString())
    ##posList = []
    ##h = 0
    ##for i in hits:
        ##posList.append(i.start())
        ##h+=1
    ##hitCount[name] = [h,posList]
    
# should find overlapping patterns too
hitCount = {}
for name in seqsDict:
    result = overlapRegEx(motifCombo,seqsDict[name].toString(),6)
    hitCount[name] = eval(result)

    x=1

hitCount_keys = hitCount.keys()
hitCount_keys.sort()

#for k in hitCount_keys:
    #print '%s\t%s' % (k,str(hitCount[k]))
    