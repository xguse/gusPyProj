import sys
import os
import weblogolib as weblogo
#from reportlab.lib import colors
#from reportlab.lib.units import cm
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Seq import UnknownSeq
#from Bio.Graphics import GenomeDiagram
#from reportlab.lib.units import cm
from gusPyCode.defs.parseSCOPE import parseScopeXML
import cPickle
from gusPyCode.defs.JamesDefs import Bag
from math import floor

def mkdirp(directory):
    if not os.path.isdir(directory):
        os.makedirs(directory)
def testForGenes(geneList, motifKey):
    geneList = set(geneList)
    hitGenes = set(motifKey[2])
    if len(geneList) == len(geneList.intersection(hitGenes)):
        return True
    else:
        return False

    

usage = '\n\npython %s outDir "geneName geneName ..." scopeOut.xml [scopeOut.xml]' % (sys.argv[0].split('/')[-1])

if len(sys.argv) < 3:
    print usage
    exit(0)


    
# 1. Open all files into own scopeOut dicts
inDict = {}
for path in sys.argv[4:]:
    inDict[path] = parseScopeXML(path)
    
# Parse GeneNames
geneNames = sys.argv[2].split()

# Print list of motifs and File names that hit genes listed:
hits = {}
for f in inDict:
    tmp = {}
    for motif in inDict[f]:
        if testForGenes(geneNames,motif):
            hits[(f,motif)] = inDict[f][motif]

hitKeys = hits.keys()
hitKeys.sort(key=lambda x: (x[0],x[1]))

if hitKeys:
    for k in hitKeys:
        if float(hits[k].sigvalue) < 1:
            del hits[k] # Remove from hits if SigVal less than 1.
        else:
            fName = k[0]
            if '/' in fName:
                fName = fName.split('/')[-1]
            print '%s\t%s\t%s\t%s\t%s\t%s' % (fName,
                                              hits[k].sigvalue,
                                              hits[k].sigValRank,
                                              hits[k].consensus,
                                              floor(len(geneNames)/float(len(hits[k].genes))*100),
                                              hits[k].algorithm)
        
    
else:
    print 'No motifs matched all of your genes. :('