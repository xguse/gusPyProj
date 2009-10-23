import cPickle
import miRNA_targeting as miRT

print '\n\n'

outFile = '/Users/biggus/Documents/James/Data/Tu_miRNA/ResultsMiRNA_Targeting/2009_10_09.seedMatches.100psCtrls.storeEvents.purged.pkl'
outFile = open(outFile,'w')
pklPath = '/Users/biggus/Documents/James/Data/Tu_miRNA/ResultsMiRNA_Targeting/2009_10_09.seedMatches.100psCtrls.storeEvents.pkl'
data    = cPickle.load(open(pklPath,'rU'))

miRs = data['miR_matches']
print 'Purging:'
for obj in miRs:
    print '\t %s' % (obj) 
    miRs[obj].purgeData()

print 'Picking purged miRNAs...'
cPickle.dump(miRs,outFile, protocol=2)

print 'Done.'


    
