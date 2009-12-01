import sys

"""Script to estimate probability of randomly choosing a misaligned nucleotide from all sequenced nucleotides
from output of yifei's SNP_Model script.  Ideally quality filtering was used to produce this script's input file.
"""

usage = 'USAGE: python estProb.py inFile'
if len(sys.argv) != 2:
    print usage
    exit(1)

totNucs = 0
varNucs = 0
    
for line in open(sys.argv[1],'rU'):
    l = line.strip('\n').split('\t')
    totNucs+=int(l[1])
    varNucs+=int(l[2])

prob = float(varNucs)/totNucs
print prob
out = open('%s.prob' % (sys.argv[1]),'w')
out.write(str(prob))