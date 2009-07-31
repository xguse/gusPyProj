from JamesDefs import revComp
from sets import Set

KmerCounts = {'AA':3,'TT':3,'GC':3,'CG':3,'CC':3,'TA':3,'AT':3}
allKmers = KmerCounts.keys()

#nrKmers = []

# make a set of sets that represent each revcomp pair Sets do not allow multiple identicle
# entries so you will end up with a single set for each rvcmp pair.  Notice that for motifs
# like 'AT' the rvcmp IS 'AT' so there is only 'AT' in the set not set([AT,AT])''. But for 
# most you should have set([fwd,rvcmp]) in rvCmpPairs
rvCmpPairSets = Set([])
for i in allKmers:
    rc = Set([i,revComp(i)])
    rvCmpPairSets.add(rc)
    
# now convert the sets back into lists and force them to become strings so we can use them for a dict key
rvCmpPairDict = {}
for i in rvCmpPairSets:
    rvCmpPairDict[(str(list(i)))] = 0
    
# iterate through the original dict and ask if the original motif key is in the list that makes up the
# rvcmp pair key.  This uses "eval" which will magicly convert the key into a list again since we created
# the string like this str([list]).  We then add the value to whichever key in rvCmpPairDict that the key
# in KmerCounts matches.
for key in KmerCounts:
    for pairKey in rvCmpPairDict:
        if key in eval(pairKey):
            rvCmpPairDict[pairKey] = rvCmpPairDict[pairKey] + KmerCounts[key] # adding kmercount to whatever value is in rvCmpPairDict[pairKey]

for each in rvCmpPairDict:
    print 'motif:%s\tcount:%s' % (each,rvCmpPairDict[each])
    
None
