import sys

assert len(sys.argv[1:]) == 2,\
       '\n\nUsage: %s inFile outFile' % (sys.argv[0].split('/')[-1])

inFile = open(sys.argv[1],'rU')
outFile = open(sys.argv[2],'w')


for line in inFile:
    data = line.strip().split('\t')
    oLine = '%s\n' % ('\t'.join([data[11],str(data[12]),str(int(data[12])+39),'>'+':'.join([data[0],data[2],data[3],data[4],data[5]]),'0','+']))
    outFile.write(oLine)