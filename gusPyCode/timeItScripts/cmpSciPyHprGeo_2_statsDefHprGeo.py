import timeit
from time import time
from gusPyCode.defs.statsDefs import hypergeoP
from scipy.stats import hypergeom

Out = '/Users/biggus/Documents/James/Data/PythonTimeTrials/20091028_scipyVsGusHprGeo_10000.txt'


runs = 100




sciTimeIt =\
"""sum([hypergeom.pmf(x, 10000, 100, 100) for x in range(3,100)])
"""
gusTimeIt =\
"""sum([hypergeoP(1000, x, 9000, 100) for x in range(3,100)])
"""


sciTimer = timeit.Timer(sciTimeIt,"from __main__ import hypergeom")
gusTimer = timeit.Timer(gusTimeIt,"from __main__ import hypergeoP")

sciTimes = sciTimer.repeat(repeat=runs,number=1)
gusTimes = gusTimer.repeat(repeat=runs,number=1)

Out = open(Out, 'w')


Out.write('scipy\tgus\n')


for i in range(runs):
    Out.write('%.6f\t%.6f\n' % (sciTimes[i],gusTimes[i]))
    
print "\n\nIM DONE!!!  "
        
        