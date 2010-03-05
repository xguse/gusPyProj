"""Create a file containing coords of rDNA exons as discribed by a GFF file."""

import csv
from gusPyCode.defs.mosqData.AaCntgCvrt import supContigConvert

gffFile = '/Users/biggus/Documents/James/Data/genomes/AaegL1/aaegypti.BASEFEATURES-AaegL1.2.gff3'
oFile   = '/Users/biggus/sandbox/testReadSubtraction/rRNAcoords.AaegL1.2.gff3.txt'
oFile   = open(oFile,'w')

reader  = csv.reader(open(gffFile, 'rU'),delimiter='\t')

for row in reader:
    if row[0].startswith('#'): 
        continue
    else:
        if row[2] == 'rRNA':
            supCntg = supContigConvert[row[0].split('|')[-1]]
            start   = row[3]
            stop    = row[4]
            ID      = row[8].split(';')[0].split('|')[-1]
            oFile.write('%s\n' % ('\t'.join([supCntg,start,stop,ID])))
oFile.flush()
oFile.close()
