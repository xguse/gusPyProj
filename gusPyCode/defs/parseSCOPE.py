import pyRXP
import xmlutils
from collections import namedtuple
import pprint
import os

#os.chdir('/Users/biggus/Documents/James/Collaborations/Campbell/scopeOut/')
# Basic Idea:
# Parse scopeXMLoutput into motif objects filterable by rank, algorithm, or genes conating the motif,
# and capable of outputing a SeqFeature(FeatureLocation()) passible to GenomeDiagram.

ltr2wrd = {'A':'adenine',
           'C':'cytosine',
           'G':'guanine',
           'T':'thymine'}


class SCOPEmotifResult(object):
    """Represents the information attributed to a single SCOPE result motif.
    
    """
    def __init__(self,xmlTagWrap_motif,rank=None):
        """Fill in soon.
        """
        xmlMotif = xmlTagWrap_motif # rename for simplicity
        
        self.rank = rank
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
        
    def _getFreq(col,base):
        """get fractional weight for a base at given column."""
        counts = [float(self.pwm['A'][col]),
                  float(self.pwm['C'][col]),
                  float(self.pwm['G'][col]),
                  float(self.pwm['T'][col])]
        tot = sum(counts)
        return float(self.pwm[base][col])/tot
    
    def toPossumMotif(self):
        """Return String representing the PossumSearch form of this motif.
        """
        
        # ++ initialize possum motif string ++
        if type(self.pwm['A'][0]) == type(1):
            pType = 'INT'
        elif type(self.pwm['A'][0]) == type(1.0):
            pType = 'FLOAT'
            
        pMtf  = 'BEGIN %s\n' % (pType)
        pMtf += 'ID %s\n'    % (self.consensus)
        pMtf += 'AC %s_%.5g\n' % (self.consensus,float(self.sigvalue))
        pMtf += 'DE %s\n'    % (self.consensus)
        pMtf += 'AL ACGT\n'
        pMtf += 'LE %s\n'    % (len(self.pwm['A']))
        
        # ++ create and add each column's data ++
        for i in range(len(self.pwm['A'])):
            pMtf += 'MA '
            for base in sorted(self.pwm.keys()):
                pMtf += '%.6g ' % (self.pwm[base][i])
            pMtf += '\n'
            
        
        pMtf += 'END\n'
        return pMtf
        
        
    def toXMSmotif(self):
        """Return string representing the XMS motif element for
        this motif.
        """
        
        # ++ initialize xms motif string ++
        xMtf = '<motif>\n\t<name>%s_%.3f</name>\n\t\t<weightmatrix alphabet="DNA" columns="%s">\n' %\
             (self.consensus,float(self.sigvalue),len(self.pwm['A']))
        # ++ create and add each column's data ++
                        
        for i in range(len(self.pwm['A'])):
            xMtf += '\t\t\t<column pos="%s">\n' % (i)
            for base in sorted(self.pwm.keys()):
                freq = self._getFreq(i,base)
                xMtf += '\t\t\t\t<weight symbol="%s">%.7f</weight>\n' % \
                     (ltr2wrd[base],freq)
            xMtf += "\t\t\t</column>\n"
        xMtf += '\t\t</weightmatrix>\n'
        xMtf += '\t\t<threshold>0</threshold>\n' # for now we'll just use null threshold
        xMtf += '\t\t<prop>\n\t\t\t<key>consensus</key>\n\t\t\t<value>%s</value></prop>\n' % \
             (self.consensus)
        xMtf += '\t\t<prop>\n\t\t\t<key>sigvalue</key>\n\t\t\t<value>%s</value></prop>\n' % \
             (self.sigvalue)
        xMtf += '\t\t<prop>\n\t\t\t<key>rank</key>\n\t\t\t<value>%s</value></prop>\n' % \
             (self.rank)
        xMtf += '\t\t<prop>\n\t\t\t<key>algorithm</key>\n\t\t\t<value>%s</value></prop>\n' % \
             (self.algorithm)
        xMtf += '</motif>\n'
        
        return xMtf
        
        
    def toTAMO(self,bkFreqs=None):
        """Return TAMO motif Object from self.pwm
        """
        
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
    
def toXMSfile(mtfDict,outPath):
    """Take a scopeMotif dict and write out the motifs in XMS format."""
    dKeys = sorted(mtfDict.keys(), key=lambda x: float(mtfDict[x].sigvalue))
    dKeys.reverse()
    oFile = open(outPath, 'w')
    oFile.write('<motifset>\n')
    for k in dKeys:
        oFile.write('%s' % (mtfDict[k].toXMSmotif()))
    oFile.write('</motifset>\n')
    oFile.close()
        
def toPossumFile(mtfDict,outPath):
    """Take a scopeMotif dict and write out the motifs in Possum format."""
    dKeys = sorted(mtfDict.keys(), key=lambda x: float(mtfDict[x].sigvalue))
    dKeys.reverse()
    oFile = open(outPath, 'w')
    oFile.write('BEGIN GROUP\n')
    for k in dKeys:
        oFile.write('%s' % (mtfDict[k].toPossumMotif()))
    oFile.write('END')
    oFile.close()

None