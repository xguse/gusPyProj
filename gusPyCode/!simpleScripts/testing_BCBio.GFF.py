"""Test BCBio.GFF
"""

from BCBio.GFF import GFFParser




gff = open('/Users/biggus/Documents/James/Data/genomes/AaegL1/aaegypti.BASEFEATURES-AaegL1.2.10000Lines.gff3')

parser = GFFParser()

#limInfo = dict(gff_type=['gene'])

for x in parser.parse(gff):
    f= x.features
    None