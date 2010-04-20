from pprint import pprint
from gusPyCode.defs.JamesDefs import pgrep

filePath = '/Users/biggus/Documents/James/Data/GeneOntology/AaGO_bp_ensembl.txt'
pattern  = 'GO:'
options  = '-c'

out = pgrep(filePath,pattern,options)
pprint(out[0:10])

None
