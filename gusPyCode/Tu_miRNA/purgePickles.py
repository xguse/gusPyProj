import cPickle
from gusPyCode.defs import miRNA_targeting as miRT

print '\n\n'

outFile = '/home/dunnw/data/tempPush/results/2009_10_26/2009_10_26.AGAP.seedMatches.100psCtrls.storeEvents.purged.pkl'
outFile = open(outFile,'w')
pklPath = '/home/dunnw/data/tempPush/results/2009_10_26/2009_10_26.AGAP.seedMatches.100psCtrls.storeEvents.pkl'
data    = cPickle.load(open(pklPath,'rU'))

miRs = data['miR_matches']
print 'Purging:'
for obj in miRs:
    print '\t %s' % (obj) 
    miRs[obj].purgeData()

print 'Picking purged miRNAs...'
cPickle.dump(miRs,outFile, protocol=2)

print 'Done.'


    
