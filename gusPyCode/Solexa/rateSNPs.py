"""Uses estimated misalign Prob to assign binomial p-Vals to each genomic position condidered and assigns 
Benjamini-Hochberg FDR corrections.
"""

from gusPyCode.defs.statsDefs import benjHochFDR
from scipy.stats import binom 
import sys

usage = '\n\n\nUSAGE: python rateSNPs.py inFile prob'
if len(sys.argv) != 3:
    print usage
    exit(1)

inFile = sys.argv[1]
prob   = float(sys.argv[2])
    
varPos = []

# Calculate and record >= x pVals
for line in open(inFile,'rU'):
    l = [int(x) for x in line.strip('\n').split('\t')]
    if l[2] > 0:
        # I think that we want cumulative p-val for x or GREATER mismatches
        #   so we use binom.cdf(x-1,n,prob) <-- need to confirm this.  
        #   *** Harsha suggests x or LESS.  I am using that untill I can ask XX.
        #       *** Tried it harshas way and things do NOT look right: 546403	29	29	0	1.0	1.0
        cumP = 1-binom.cdf(l[2]-1,l[1],prob)
        varPos.append(l+[cumP])

# Calculate BH adjusted q-vals
varPos = benjHochFDR(varPos,pValColumn=4)

fOut = open('%s.pr%.5f.qVals' % (inFile,prob), 'w')
for item in varPos:
    fOut.write('%s\n' % ('\t'.join([str(x) for x in item])))
        
        