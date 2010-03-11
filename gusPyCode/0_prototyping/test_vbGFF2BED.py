from gusPyCode.defs.bioDefs import vbGFF2BED

gff = '/Users/biggus/Documents/James/Data/genomes/AaegL1/aaegypti.BASEFEATURES-AaegL1.2.1000Lines.gff3'
bed = '/Users/biggus/Documents/James/Data/genomes/AaegL1/testGFF2BED.txt'


vbGFF2BED(gff,bed,filterOnType='exon')