print 'loading modules...'
import motility
from gusPyCode.gClasses.gSeqClasses import DNAseq

print 'loading seq file...'
seqFile    = map(lambda line: line.strip(), open('/Users/biggus/Documents/James/Data/AedesCorePromoterWork/output/TATA_Inr/Aa_Ensbl49_AaegL1.1.plus50minus100._TataInr_seqSlices_.fas','rU').readlines())
outFile    = '/Users/biggus/Documents/James/Data/AedesCorePromoterWork/output/TATA_Inr/Aa_Ensbl49_AaegL1.1.plus50minus100._TataInr_seqSlices_.locations.txt'
motif2 = motility.IUPAC('RGWYV') 
winTo2ndMotif = (23,27) # <---- This is based from END of first motif!!!
seqList = zipped = zip(seqFile[:-1:2], seqFile[1::2])
seqDict = {}

print 'populating seqDict...'
for item in seqList:
    name      = item[0].replace('|Inr(..DPE)','').replace('>','')
    header    = item[0]
    seq       = item[1]
    coverage  = item[0].split('|')[-1].split('-')
    coverage  = (int(coverage[0]),int(coverage[1]))
    motif1Loc = coverage[0]-6
    
    seqDict[name]           = DNAseq(seq,name,header)
    seqDict[name].coverage  = coverage
    seqDict[name].motif1Loc = motif1Loc


def checkSlice(seqObj, motif2, winTo2ndMotif):
    motif2Loc = None
    mHits = motif2.find(seqObj.toString())
    
    for each in mHits:
        if each[2] == -1:
            continue
        elif each[0] >= winTo2ndMotif[0] and each[0]<= winTo2ndMotif[1]:
            motif2Loc = each[0]
    if motif2Loc:
        return '%s\t%s\t%s\t%s\n' % (seqObj.name,seqObj.motif1Loc,motif2Loc+seqObj.coverage[0],motif2Loc)
    else:
        return None
    
seqDict_sortedKeys = seqDict.keys()
seqDict_sortedKeys.sort()

resultList = []
print 'checking slices...'
for name in seqDict_sortedKeys:
    result = checkSlice(seqDict[name], motif2, winTo2ndMotif)
    
    if result:
        resultList.append(result)
    else:
        continue
    
outFile = open(outFile,'w')
print 'writing to outFile...'
outFile.writelines(resultList)

print 'Done.'
    

    

