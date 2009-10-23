import motility
from gusPyCode.gClasses.gSeqClasses import DNAseq


class HitList:
    def __init__(self,seqHitsDict):
        self.data = seqHitsDict
        self.sortedKeys = seqHitsDict.keys()
        self.sortedKeys.sort()
        self.seqsWithHit = {}
        self.countOfSeqsWithHit = None
        
    def countSeqWithHit(self):
        count = 0
        for key in self.sortedKeys:
            for hit in self.data[key]:
                if hit[2] == 1:
                    count+=1
                    self.seqsWithHit[key] = self.data[key]
                    ##print self.data[key]
                    break
                
        self.countOfSeqsWithHit = count
print 'loading file...'
seqFile    = map(lambda line: line.strip(), open('/Users/biggus/Documents/James/Data/AedesCorePromoterWork/sourceCoords/Aa_Ensbl49_AaegL1.1.Plus50Minus100coords.geneStrand.fas','rU').readlines())
outFile    = '/Users/biggus/Documents/James/Data/AedesCorePromoterWork/output/Aa_Ensbl49_AaegL1.1.plus50minus100._InrDPE_Hits_2.fas'
outFile2   = '/Users/biggus/Documents/James/Data/AedesCorePromoterWork/output/Aa_Ensbl49_AaegL1.1.plus50minus100._InrDPE_seqSlices_2.fas'
motif1 = motility.IUPAC('TCAKTY')
motif2 = motility.IUPAC('RGWYV')
winTo2ndMotif = 36
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
print 'searching seqDict...'
for each in seqNames:
    seqHits[each] = motif1.find(seqDict[each].toString())

### write to file in an easily python readable format
##outFile = open(outFile, 'w')

##for every in seqNames:
    ##outFile.write('%s:%s\n' % (every,seqHits[every]))
print 'creating hitList1...'
hitList1  = HitList(seqHits)

hitList1.countSeqWithHit()
##print "I found %s seqs with motif." % (hitList.countOfSeqsWithHit)
##print 'len of seqsWithHit dict = %s' % (len(hitList.seqsWithHit))

# create new seqDict for seqs with motif1
hitSeqDict = {}
motif1SeqsNames = hitList1.seqsWithHit.keys()

print 'transfering hits from 1st motif to hitsDict...'
for each in motif1SeqsNames:
    hitSeqDict[each] = seqDict[each]
print 'hitSeqDict len = %s...' % (len(hitSeqDict))

hitSeqDict_keys  = hitSeqDict.keys()
doubleHitSeqDict = {}
seqSliceDict = {}

for key in hitSeqDict_keys:
    for site in hitList1.seqsWithHit[key]:
        # creat slice and attach info to it
        seqSlice = seqDict[key].slice(int(site[1]),int(site[1])+winTo2ndMotif)
        seqSlice.name = '%s|%s(..%s)|%s-%s' % (seqDict[key].name,motif1.motif,motif2.motif,site[1],int(site[1])+winTo2ndMotif)
        ##print seqSlice.name
        # ommit seqs that are shorter than motif to be searched
        if len(seqSlice) < len(motif2):
            continue
        else:
            
            motif2hits  = motif2.find(seqSlice.toString())
        
        for hit in motif2hits:
            if hit[2] == 1:                
                ##print "%s:%s" % (key,hit)
                # on success drop SeqObj into doubleHitSeqDict
                doubleHitSeqDict[key] = seqDict[key]
                # drop actual seq into seqSliceDict
                seqSliceDict[seqSlice.name] = seqSlice
                print '%s->%s' % (seqSlice.name,len(seqSlice))
                #break
print 'doubleHitSeqDict len = %s\nseqSliceDict = %s' % (len(doubleHitSeqDict), len(seqSliceDict))

print 'writing double hits to fasta file...'
outFile = open(outFile,'w')
# sort keys
doubleHitSeqDict_keys = doubleHitSeqDict.keys()
doubleHitSeqDict_keys.sort()
for name in doubleHitSeqDict_keys:
    outFile.write(doubleHitSeqDict[name].toFasta())

print 'writing seqSlices to fasta file...'
outFile2 = open(outFile2,'w')
# sort keys
seqSliceDict_keys = seqSliceDict.keys()
seqSliceDict_keys.sort()
for name in seqSliceDict_keys:
    outFile2.write(seqSliceDict[name].toFasta())


print 'Done.'