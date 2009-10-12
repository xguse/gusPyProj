import sys

"""
Extracts the unique ID and real seq from the file type: s_N_sequence.txt.
expl:
@SOLEXA2:7:1:1:1999#0/1
NAAAAGCGAGATCGCGATCGAAGGCGGTNNNNNNNNNNNN
+SOLEXA2:7:1:1:1999#0/1
DLUQUUVUUUUVWUVUBBBBBBBBBBBBBBBBBBBBBBBB
"""
args = sys.argv[:]

usage = "USAGE: %s inFile" % (args[0])

assert len(args) == 2, usage

inFile = open(args[1], 'rU')

##inFile = '/Users/biggus/Downloads/s_7_sequence.500Lines.txt'
##inFile = open(inFile, 'rU')

# Read in file one line at a atime and check for a '@' at the start
# If so: ignor info after @ and read in next line (this is the seq)

seqDict = {}

while 1:
    line = inFile.readline()
    if not line: break
    if line.startswith('@'):
        seq = inFile.readline().strip('\n')
        if seq in seqDict:
            seqDict[seq] += 1
        else:
            seqDict[seq] = 1

count = 0
for seq in seqDict:
    count+=1
    print '>%08.f_%s\n%s' % (count,seqDict[seq], seq)
    
