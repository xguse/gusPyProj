import sys


usage = "USAGE: python %s <inFile> <outFile>" % (sys.argv[0].split('/')[-1])

assert len(sys.argv) == 3, usage

inFile  = sys.argv[1]
outFile = sys.argv[2]
inFile  = open(inFile,'r')

counts = []

while 1:
    line = inFile.readline()
    if not line: break
    if line.startswith('>'):
        counts.append(line.strip().split('_')[-1]+'\n')
        
open(outFile, 'w').writelines(counts)

