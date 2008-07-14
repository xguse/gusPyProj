import JamesDefs

from Bio import Seq
import os
#motifList = ['[ta]tc[ctg]ttt', 'tt[ag]g[tc]at]']
#motifsInAll = countMotifInAll(motifStr, seqDict)

d = {'a':a = Seq('atccttt'), 'b':Seq('ttcgttt'), 'c':Seq('aaaggata'), 'd': Seq('AAAGGATA')}


for k in d:
    print d[k]

print 'yay'