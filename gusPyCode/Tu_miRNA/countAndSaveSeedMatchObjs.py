import JamesDefs
import miRNA_targeting as miTrgt
import sys
import cPickle
from time import time

#test

outFile = '/home/dunnw/data/tempPush/results/2009_10_26/2009_10_26.AGAP.seedMatches.100psCtrls.storeEvents.pkl'
     #'/Users/biggus/Documents/James/Data/Tu_miRNA/ResultsMiRNA_Targeting/DATE.seedMatches.NUMandTYPEofCTRLS.pkl'

initCtrlsWith = 'proSeed' # 'proSeed' or 'matchVersion'     
ctrlIter      = 100 
refGenome     = 'AGAP'


print '-- -- -- -- --\n\nLoading files...'
orthoPath  = '/home/dunnw/data/tempPush/3wayOrthos.combineOrthosJamesDefs.out.filtered.noInfers.txt'
seqPaths   = ['/home/dunnw/data/tempPush/Aa_500afterCoding.usuable.stpCdn.fas',
              '/home/dunnw/data/tempPush/Ag_500afterCoding.usuable.stpCdn.fas',
              '/home/dunnw/data/tempPush/Cq_500afterCoding.newCoords.usuable.stpCdn.fas']          
miRNA_Path = '/home/dunnw/data/tempPush/mature.aga.fa'

          
seqDict        = miTrgt.loadSeqs(seqPaths)
orthoRelations = miTrgt.loadOrthos(orthoPath, seqDict)
miRNAs         = miTrgt.loadMiRNAs(miRNA_Path)

# randGroup = ['aga-miR-12','aga-miR-263b']
# d = {}
# for i in randGroup:
#     d[i] = miRNAs[i]
# miRNAs = d

print 'Getting ortho seqs...'
orthoSeqs  = miTrgt.filterOrthoSeqs(seqDict,orthoRelations)



seenSeeds = set()
miR_matches = {}

saveObj = {'orthoRelations':orthoRelations,
           'seenSeeds':seenSeeds,
           'miR_matches':miR_matches}


print 'Initializing matchVersions...'

for m in miRNAs:
    seed = miTrgt.miRNA(miRNAs[m],seenSeeds,orthoRelations=orthoRelations,name=m)
    miR_matches[seed.name]= seed

print 'Initializing ctrls with %s...' % (initCtrlsWith)#'Initializing ctrls from MatchVersion...'
if initCtrlsWith == 'matchVersion':
    for m in miR_matches:
        miR_matches[m].buildCtrlsFromMatchVers(seenSeeds,ctrlIter) #miR_matches[m].buildCtrlsFromProSeed(seenSeeds,100)
elif initCtrlsWith == 'proSeed':
    for m in miR_matches:
        miR_matches[m].buildCtrlsFromProSeed(seenSeeds,ctrlIter) #miR_matches[m].buildCtrlsFromProSeed(seenSeeds,100)

t1 = time()
print 'Tallying hits...'
for m in miR_matches:
    miR_matches[m].tallyHits(orthoSeqs)
    print '\t'+m

if refGenome:
    print 'Counting hits in orthos for %s...' % (refGenome)
else:
    print 'Counting hits in orthos...'
for m in miR_matches:
    if refGenome:
        miR_matches[m].countHitsInOrthos4(genomeToken=refGenome,returnGenes=True)
    else:
        miR_matches[m].countHitsInOrthos(returnGenes=True)
t2=time()
print 'Tallying and counting took %.2f min.' % ((t2-t1)/60.0)

print 'Pickling saveObj...'
outFile = open(outFile, 'w')
cPickle.dump(saveObj,outFile, protocol=2)

print 'Done.'
