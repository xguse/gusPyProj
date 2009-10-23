import motility
from gSeqClasses import DNAseq


class HitList:
    def __init__(self,seqHitsDict):
        self.data = seqHitsDict
        self.sortedKeys = seqHitsDict.keys()
        self.sortedKeys.sort()
        self.seqsWithHit = None
        
    def countSeqWithHit(self):
        count = 0
        for key in self.sortedKeys:
            for hit in self.data[key]:
                if hit[2] == 1:
                    count+=1
                    break
                
        self.seqsWithHit = count

seqFile    = map(lambda line: line.strip(), open('/Users/biggus/Documents/James/Data/AedesCorePromoterWork/sourceCoords/Aa_Ensbl49_AaegL1.1.Plus50Minus100coords.geneStrand.fas','rU').readlines())
outFile    = '/Users/biggus/Documents/James/Data/AedesCorePromoterWork/Aa_Ensbl49_AaegL1.1.plus50minus100._InrHits_.txt'
iupacMotif = motility.IUPAC('TCAKTY')
pfmMotif   = motility.PWM([[196,230,260,456],[54,83,94,911],[20,963,18,141],[1052,13,16,61],[24,50,692,376],[27,21,3,1091],[43,440,60,599],[206,216,407,313],[251,278,294,319],[285,285,197,375],[255,251,228,408],[259,207,237,439],])

seqList = zipped = zip(seqFile[:-1:2], seqFile[1::2])

seqDict = {}

# populate seqDict
for item in seqList:
    name   = item[0].split(' ')[0].split(':')[0].replace('>','')
    header = item[0]
    seq    = item[1]
    
    seqDict[name] = DNAseq(seq,name,header)
    

seqHits  = {}
seqNames = seqDict.keys()
seqNames.sort()

# populate Hits Dict
for each in seqNames:
    #seqHits[each] = iupacMotif.find(seqDict[each].toString(),1)
    seqHits[each] = pfmMotif.find(seqDict[each].toString(),6600)
### write to file in an easily python readable format
##outFile = open(outFile, 'w')

##for every in seqNames:
    ##outFile.write('%s:%s\n' % (every,seqHits[every]))

hitList  = HitList(seqHits)

hitList.countSeqWithHit()
print "I found %s seqs with motif." % (hitList.seqsWithHit)

print 'Done.'