import cPickle
import miRNA_targeting as miRT

print '\n\n'


pklPath = '/Users/biggus/Documents/James/Data/Tu_miRNA/ResultsMiRNA_Targeting/2009_08_31.seedMatches.100ctrls.pkl'
data    = cPickle.load(open(pklPath,'rU'))


print 'Im Done.'