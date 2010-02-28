from gusPyCode.defs.bioDefs import ParseFastQ

fqFile = '/Users/biggus/sandbox/testParseFastQ/fastq.txt'

p = ParseFastQ(fqFile)

while 1:
    elem = p.getNext()
    if not elem:
        break
    else:
        print elem