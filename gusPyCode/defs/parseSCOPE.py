import pyRXP
import xmlutils
from collections import namedtuple
import pprint
import os

os.chdir('/Users/biggus/Documents/James/Collaborations/Campbell/scopeOut/')
# Basic Idea:
# Parse scopeXMLoutput into motif objects filterable by rank, algorithm, or genes conating the motif,
# and capable of outputing a SeqFeature(FeatureLocation()) passible to GenomeDiagram.



class SCOPEmotifResult(object):
    """Represents the information attributed to a single SCOPE result motif.
    
    """
    def __init__(self,xmlTagWrap_motif,rank):
        """Fill in soon.
        """
        xmlMotif = xmlTagWrap_motif # rename for simplicity
        
        self.sigValRank = rank
        self.consensus = xmlMotif.sequence
        self.sigvalue = xmlMotif.sigvalue
        self.algorithm = xmlMotif.algorithm
        self.genes = tuple(sorted(list(set([x.gene for x in xmlMotif.instances])))) # (sorted tuple)
        self.instances = () # tuple of namedtuples keys=(sequence, strand<-/+>,start[negVersion],end[negVersion],geneName))
        self.pwm = {} # key(nulceotide),value(tuple of freqs at each position)
        
        # Build instances namedtuples:
        Instance = namedtuple('instance','sequence strand begin end gene')
        instList = [Instance(x.sequence, x.strand, x.begin, x.end, x.gene) for x in xmlMotif.instances]
        instList = tuple(instList)
        self.instances = instList
        
        # Build PWM dict:
        for base in xmlMotif.pwm:
            self.pwm[base.name.upper()] = tuple([int(x[2][0]) for x in base._children])
        
    
    def toTAMO(self,bkFreqs=None):
        """Return TAMO motif Object from self.pwm
        """
        pass
    def toTRANSFACstring(self):
        """Returns printable motif representation in TRANSFAC format sutible for use in MAST.
        """
        pass
    def toGDiagramFeature(self, geneName):
        """Return a list of all seqFeature(FeatureLocation()) objs for a given geneName.
        """
        
        pass

def parseScopeXML(pathToXML):
    """For a scopeXML file, return a dict of SCOPE result objects with key defined as a tuple of
    the motif's sigValue rank<int>,the algorithm<string>, and a <tuple> of the genes it was found in:
    (rank,algorithm,(gene,gene,...,gene))
    This will allow sorting and filtering on the key for these values.
    """
    xmlString = ''.join(map(lambda l: l.strip(), open(pathToXML, 'rU')))
    tree = pyRXP.Parser().parse(xmlString)
    tagWrap = xmlutils.TagWrapper(tree)
    
    # collect which genes were used and how long they were:
    analyzedGeneInfo = {} # key(geneName), value(length of promoter considered)
    for g in tagWrap.genes:
        analyzedGeneInfo[g.name] = g.length
        
    rDict= {}
    for i in range(len(tagWrap.motifs)):
        genes = tuple(sorted(list(set([x.gene for x in tagWrap.motifs[i].instances]))))
        rank  = i+1
        rDict[(rank,tagWrap.motifs[i].algorithm,genes)] = SCOPEmotifResult(tagWrap.motifs[i],rank)
    
    return rDict
    
    
None