import timeit
from gmpy import comb
from gusPyCode.defs.doug_hypergeometric import binc
from gusPyCode.defs.statsDefs import hypergeoP,hypergeoP_Test

n=40
i=25
m=60
N=30

testRuns = 100000

methods = {}

#methods['sdBinom']  = timeit.Timer('binom(%s,%s)' % (n,k), 'from gusPyCode.defs.statsDefs import binom')
#methods['sdChoose']  = timeit.Timer('choose(%s,%s)' % (n,k), 'from gusPyCode.defs.statsDefs import choose')
methods['hG']   = timeit.Timer('hypergeoP(%s,%s,%s,%s)'% (n,i,m,N), 'from gusPyCode.defs.statsDefs import hypergeoP')
methods['hG_T'] = timeit.Timer('hypergeoP_Test(%s,%s,%s,%s)'% (n,i,m,N), 'from gusPyCode.defs.statsDefs import hypergeoP_Test')

for k in sorted(methods.keys()):
    print 'Method(%s) took %.5f seconds.' % (k,methods[k].timeit(testRuns))

 

