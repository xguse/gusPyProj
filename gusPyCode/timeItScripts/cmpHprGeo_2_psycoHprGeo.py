
import timeit
from time import time
from gusPyCode.defs.statsDefs import hypergeoP
#from scipy.stats import hypergeom

Out = '/Users/biggus/Documents/James/Data/PythonTimeTrials/20091031_gusHprGeoVsPsyco_100.txt'


runs = 100




stdTimeIt =\
"""sum([hypergeoP(1000, x, 9000, 100) for x in range(3,100)])
"""
#"""sum([hypergeom.pmf(x, 10000, 100, 100) for x in range(3,100)])
#"""

psyTimeIt =\
"""sum([hypergeoP(1000, x, 9000, 100) for x in range(3,100)])
"""
proTimeIt =\
"""sum([hypergeoP(1000, x, 9000, 100) for x in range(3,100)])
"""

stdTimer = timeit.Timer(stdTimeIt,"from __main__ import hypergeoP")
psyTimer = timeit.Timer(psyTimeIt,"from __main__ import hypergeoP\nimport psyco\npsyco.full()")
proTimer = timeit.Timer(proTimeIt,"from __main__ import hypergeoP\nimport psyco\npsyco.profile()")

stdTimes = stdTimer.repeat(repeat=runs,number=1)
psyTimes = psyTimer.repeat(repeat=runs,number=1)
proTimes = proTimer.repeat(repeat=runs,number=1)

Out = open(Out, 'w')


Out.write('stdpy\tpsy\tpro\n')


for i in range(runs):
    Out.write('%.6f\t%.6f\t%.6f\n' % (stdTimes[i],psyTimes[i],proTimes[i]))
    
print "\n\nIM DONE!!!  "
        
        