import cPickle
from gusPyCode.defs import miRNA_targeting as miRT

print '\n\n'


pklPath = '/Users/biggus/Documents/James/Data/Tu_miRNA/ResultsMiRNA_Targeting/test.seedMatches.100psCtrls.storeEvents.pkl'
data    = cPickle.load(open(pklPath,'rU'))


print 'Im Done.'