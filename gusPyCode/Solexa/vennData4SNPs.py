from time import time
import optparse
import sys

from gusPyCode.defs.mosqData.AaCntgCvrt import supContigConvert


def printSNPlines(exons,strainLines):
    
    for ex in exons:
        for line in strainLines:
            if line[1] == ex:
                for snp in eval(line[-1]):
                    print '\t'.join([ex]+snp)
                print ''


topX = 50


Liv = 'LSandB.minus_rRNA.varscan.snps.combined.readcounts.filtered.snpsPerBp.cov50.txt'
Rex = 'RBandS.minus_rRNA.varscan.snps.combined.readcounts.filtered.snpsPerBp.cov50.txt'
Che = 'CBandS.minus_rRNA.varscan.snps.combined.readcounts.filtered.snpsPerBp.cov50.txt'

Liv = map(lambda l: l.strip('\n').split('\t'), open(Liv,'rU'))
Rex = map(lambda l: l.strip('\n').split('\t'), open(Rex,'rU'))
Che = map(lambda l: l.strip('\n').split('\t'), open(Che,'rU'))

Che.pop(0);Liv.pop(0);Rex.pop(0)

Che.sort(key=lambda x: float(x[5]))
Liv.sort(key=lambda x: float(x[5]))
Rex.sort(key=lambda x: float(x[5]))

Che.reverse()
Liv.reverse()
Rex.reverse()


CXexns = set([x[1] for x in Che[:topX]])
LXexns = set([x[1] for x in Liv[:topX]])
RXexns = set([x[1] for x in Rex[:topX]])

CLRX = CXexns.intersection(LXexns.intersection(RXexns))
CLX  = CXexns.intersection(LXexns).difference(CLRX)
CRX  = CXexns.intersection(RXexns).difference(CLRX)
LRX  = LXexns.intersection(RXexns).difference(CLRX)

print \
'''Here are the Venn Data:
CLRX: %s
CLX: %s
CRX: %s
LRX: %s''' % (len(CLRX),
              len(CLX),
              len(CRX),
              len(LRX))

print "--=[CLRX]=--"
print "Chetumal:"
printSNPlines(list(CLRX),Che)
print "\nLiverpool:"
printSNPlines(list(CLRX),Liv)
print "\nRexD:"
printSNPlines(list(CLRX),Liv)

print "--=[CLX]=--"
print "Chetumal:"
printSNPlines(list(CLX),Che)
print "\nLiverpool:"
printSNPlines(list(CLX),Liv)

print "--=[CRX]=--"
print "Chetumal:"
printSNPlines(list(CRX),Che)
print "\nRexD:"
printSNPlines(list(CRX),Liv)

print "--=[LRX]=--"
print "Liverpool:"
printSNPlines(list(LRX),Liv)
print "\nRexD:"
printSNPlines(list(LRX),Liv)
