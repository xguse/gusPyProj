#!/usr/bin/env python

import sys

args = sys.argv[1:]

usage = """\n\nUSAGE: %s inFile outFile""" % (sys.argv[0].split('/')[-1])

assert len(args) == 2, usage

inPath  = args[0]
outPath = args[1]
inFile  = open(inPath,'rU')
outFile = open(outPath,'w')

contigs      = {}
contigCounts = {}

# Step 1: Cycle Through file collecting a list of all contigs presnt in the file and
#         cataloging how many times an exact read is encountered for the contig.
print 'Tallying data...'
for line in inFile:
    fields = line.strip('\n').split('\t')    # [Q1]
    if fields[11] in contigs:                # [Q2.a]
        contigCounts[fields[11]] += 1        # [Q2.a]
    else:
        contigCounts[fields[11]] = 1         # [Q2.a]
        contigs[fields[11]] = {}             # [Q2.a]
    
    if fields[8] in contigs[fields[11]]:     # [Q2.b]
        contigs[fields[11]][fields[8]] += 1  # [Q2.b]
    else:
        contigs[fields[11]][fields[8]] = 1   # [Q2.b]
    
        
def getNumUniqReads(contigName,contigsDict):          # [Q3]
    uniqs = []                                        # [Q3]
    for readSeq in contigsDict[contigName]:           # [Q3]
        if contigsDict[contigName][readSeq] == 1:     # [Q3]
            uniqs.append(readSeq)                     # [Q3]
    return tuple(uniqs)                               # [Q3]
        
# Step 2: write the results out to the outFile
print 'Writing to file...'        
for contigName in sorted(contigs.keys()):
    outFile.write('%s\t%s\t%s\n' \
                  % (contigName,
                     contigCounts[contigName],
                     len(getNumUniqReads(contigName,contigs))))

print 'Im done!'