import sys
import random
from scipy import mat, transpose
from TAMO.MotifTools import Motif
from TAMO.seq import Fasta
from gusPyCode.defs.statsDefs import hypergeoP

class ParseSolSorted(object):
    """Class to parse and return a single read entry from solexa x_sorted.txt file type."""
    def __init__(self,filePath):
        """Returns a line-by-line solexa x_sorted.txt parser analogous to file.readline().
        Exmpl: parser.getNext() """
        self._file = open(filePath, 'rU')
        
    def getNext(self):
        """Reads in next line, parses fields, returns fieldTuple or None (eof)."""
        line = self._file.readline()
        if line:
            return tuple(line.strip('\n').split('\t'))
        else: 
            return None
    
    def getNextReadSeq(self):
        """Calls self.getNext and returns only the readSeq."""
        line = self.getNext()
        if line:
            return line[8]
        
    def getNextReadCoords():
        """Calls self.getNext and returns only the readCoords t(str(contig),int(start),int(stop))."""
        line = self.getNext()
        if line:
            contig = line[11]
            start  = int(line[12])
            stop   = int(line[12])+int(line[14])-1 # start+len-1
            return tuple([contig,start,stop])


class ParseBowtieBed(object):
    """Class to parse and return a single read entry from bowtie_bed file type."""
    def __init__(self,filePath):
        """Returns a line-by-line bowtie_bed parser analogous to file.readline().
        Exmpl: parser.getNext() """
        self._file = open(filePath, 'rU')
        
    def getNext(self):
        """Reads in next line, parses fields, returns fieldTuple or None (eof)."""
        line = self._file.readline()
        if line:
            return tuple(line.strip('\n').split('\t'))
        else: 
            return None
    
    def getNextReadSeq(self):
        """Calls self.getNext and returns only the readSeq."""
        line = self.getNext()
        if line:
            return line[3].split('_')[-1]


class ParseFastQ(object):
    """Returns a read-by-read fastQ parser analogous to file.readline()"""
    def __init__(self,filePath):
        """Returns a read-by-read fastQ parser analogous to file.readline().
        Exmpl: parser.getNext()"""
        self._file = open(filePath, 'rU')
        self._currentLineNumber = 0
    
    def getNext(self):
        """Reads in next element, parses, and does minimal verification.
        Returns: tuple: (seqHeader,seqStr,qualHeader,qualStr)"""
        # ++++ Get Next Four Lines ++++
        elemList = []
        for i in range(4):
            line = self._file.readline()
            self._currentLineNumber += 1 ## increment file position
            if line:
                elemList.append(line.strip('\n'))
            else: return None
        
        # ++++ Check Lines For Expected Form ++++
        # -- Make sure we got 4 full lines of data --
        assert "" not in elemList,\
               "** ERROR: It looks like I encountered a premature EOF or empty line.\n\
               Please check FastQ file near line #%s (plus or minus ~4 lines) and try again**" % (self._currentLineNumber)
        # -- Make sure we are in the correct "register" --
        assert elemList[0].startswith('@'),\
               "** ERROR: The 1st line in fastq element does not start with '@'.\n\
               Please check FastQ file and try again**"
        assert elemList[2].startswith('+'),\
               "** ERROR: The 3rd line in fastq element does not start with '+'.\n\
               Please check FastQ file and try again**"
        
        # ++++ Return fatsQ data as tuple ++++
        return tuple(elemList)
    
    def getNextReadSeq(self):
        """Calls self.getNext and returns only the readSeq."""
        record = self.getNext()
        if record:
            return record[1]
        else:return None


class MASTmotif:
    """Represents MAST type motifs."""
    def __init__(self,TAMOmotif):
        self.alpha = 'ACGT'
        self.alen  = 4
        self.w     = len(TAMOmotif)
        self.mat   = []
        
        for pos in TAMOmotif.logP:
            self.mat.append([pos['A'],pos['C'],pos['G'],pos['T']])
    
    def toText(self):
        """Returns string form of the motif in MAST format."""
        
        out = 'log-odds matrix: alength= 4 w= %s\n' % (self.w)
        for pos in self.mat:
            out+='%s\t%s\t%s\t%s\n' % (pos[0],pos[1],pos[2],pos[3])
        return out


def goEnrichment(geneCluster,GOgenes,popSize):
    """
    goEnrichment(geneCluster,GOgenes,popSize):
    geneCluster = set(genesGroupedBySomeQuality)
    GOgenes     = set(genesInGOgroup)
    popSize     = int(numberOfGenesConsideredAsPopulation)
    
    Returns cumulative hypergeometric enrichment P-value.
    """
    #n = # of positives in population
    #i = # of positives in sample
    #m = # of negatives in population
    #N = sample size
    #P(x=i) = (choose(n,i)choose(m,N-i))/choose(n+m,N)
    #For more details -> http://mathworld.wolfram.com/HypergeometricDistribution.html
    n = len(GOgenes)
    i = len(GOgenes.intersection(geneCluster))
    m = popSize-len(GOgenes)
    N = len(geneCluster)
    
    return sum([hypergeoP(n,x,m,N) for x in range(i,N+1)])


def makeRandomSeq(length):
    alpha = ['A','C','G','T']
    seq   = '' 
    
    for i in range(length):
        seq = seq + random.choice(alpha)
    
    return seq

def taMotif2MopatMatrix(taMotif):
    """
    Uses ONE tamo motif to produce a probability matrix compatible with MOPAT.
    RETURNS: [oneLetterName, [list of lists]]
    """
    name = taMotif.oneletter.replace('.','n')
    matrix = []
    
    # Use logLiklihoods to generate info for motif.counts
    taMotif = Motif(taMotif.bogus_kmers())
    
    for pos in taMotif.counts:
        mopatPos = []
        for nuc in sorted(pos.keys()):
            mopatPos.append(pos[nuc])
        matrix.append(mopatPos)
            
    return [name, matrix]

def convertTAMOs2MOPATs(listOfTamoMotifs):
    """
    Converts a list of tamo motifs into a dict of mopat compatible matrices.
    RETURNS: dictOfMopatMatrices (keys<number_oneletter>, values<listOfLists>)
    """
    
    mopatDict = {}
    
    count = 0
    for m in listOfTamoMotifs:
        count+=1
        temp = taMotif2MopatMatrix(m)
        mopatDict['%s_%s' % (count,temp[0])] = temp[1]
    
    return mopatDict

#def jaspar2tamoMotif(jasparMatrix):
    
    ## Transpose the matrix
    #pfm = transpose(mat(jasparMatrix)).tolist()
    
    ## convert to tamo-like count dict
    
    
    
    #pass

def softMaskDict2HardMask(fastaDict):
    """
    Replaces(in place) all lower-case sequence letters with 'N'.
    """
    
    for i in fastaDict:
        fastaDict[i] = fastaDict[i].replace('a','N')
        fastaDict[i] = fastaDict[i].replace('c','N')
        fastaDict[i] = fastaDict[i].replace('g','N')
        fastaDict[i] = fastaDict[i].replace('t','N')
    

def geneList2FastaDict(geneList, sourceFastaPath, hardMasked=True):
    """
    Returns a Dict of requested fasta recs in form SeqName:Sequence.
    Defaults to HardMasked return seqeunces.
    """
    
    sourceDict = Fasta.load(sourceFastaPath)
    
    # make new dict of all genes both in geneList AND sourceDict
    # new dict may be shorter than geneList!!!!!!
    
    newDict = {}
    for i in geneList:
        if sourceDict[i]:
            newDict[i] = sourceDict[i]
            
    print "%s genes names given, %s found." % (len(geneList), len(newDict))
    
    if hardMasked:
        softMaskDict2HardMask(newDict)
    
    return newDict

def ifKmerInAll(kmer,dictOfSeqs, factor=1):
    kmer = Motif(kmer)
    results = []
    for seq in dictOfSeqs:
        temp = kmer.scan(dictOfSeqs[seq], factor=factor)
        if temp[0]:
            results.append(True)
        else:
            results.append(False)
        
    if False in results:
        return False
    else:
        return True

def bestIdentOverLen(Bio_Blast_Record_Blast): # full name of BioPython obj for clairity
    """Returns the greatest fraction of the query's length found to be identicle
    in any alignment objs HSPs. <possible values: 0 to 1>"""

    blastRec = Bio_Blast_Record_Blast # rename for ease of use
    
    bestIdent = 0
    for i in range(len(blastRec.alignments)):
        # use alignments[].hsps[0].match to count how many consequtive |s
        None
        
    # !!! this def is NOT finished !!! #
    assert 1==2, ' def bestIdentOverLen() is NOT finished!!!'
    

def fastaFileToBioSeqDict(pathToFastaFile, Alphabet='IUPACAmbiguousDNA', splitOn=None, indexAsName=0):
    from Bio.Alphabet.IUPAC import Alphabet
    from Bio import SeqIO
    #  Populate a Dict with Seq objs for Anopheles boundary seqs
    #  What follows directly is a klugde to get all my seqDict vals to have an IUPAC alphabet
    listOfSeqs = list(SeqIO.parse(open(pathToFastaFile, "rU"), "fasta"))
    for record in listOfSeqs :
        record.seq.alphabet = Alphabet

    dictOfSeqs = SeqIO.to_dict(listOfSeqs, key_function = lambda rec : rec.description.split(splitOn)[indexAsName])

    return dictOfSeqs

#=========================================================================
compl_iupacdict = {'A':'T',
                   'C':'G',
                   'G':'C',
                   'T':'A',
                   'M':'K',
                   'R':'Y',
                   'W':'W',
                   'S':'S',
                   'Y':'R',
                   'K':'M',
                   'V':'B',
                   'H':'D',
                   'D':'H',
                   'B':'V',
                   'X':'X',
                   'N':'N',
                   'a':'t',
                   'c':'g',
                   'g':'c',
                   't':'a',
                   'm':'k',
                   'r':'y',
                   'w':'w',
                   's':'s',
                   'y':'r',
                   'k':'m',
                   'v':'b',
                   'h':'d',
                   'd':'h',
                   'b':'v',
                   'x':'x',
                   'n':'n',
                   '-':'-'}

def compliment(motif, compl_iupacdict):
    compl_motif = ""
    for i in range(0,len(motif)):
        letter = motif[i]
        compl_motif = compl_motif + compl_iupacdict[letter]
    return compl_motif

def reverse(text):
    return text[::-1]

def revComp(seq):
    revCompSeq = reverse(compliment(seq, compl_iupacdict))
    return revCompSeq
#=========================================================================

def iupacList_2_regExList(motifList):
    i = 0
    while i < len(motifList):
        motifList[i] = [motifList[i], iupac2regex(motifList[i])]
        i += 1




def iupac2regex(motif):

    iupacdict = {'A':'A',
                 'C':'C',
                 'G':'G',
                 'T':'T',
                 'M':'[AC]',
                 'R':'[AG]',
                 'W':'[AT]',
                 'S':'[CG]',
                 'Y':'[CT]',
                 'K':'[GT]',
                 'V':'[ACG]',
                 'H':'[ACT]',
                 'D':'[AGT]',
                 'B':'[CGT]',
                 'X':'[ACGT]',
                 'N':'[ACGT]'}

    transl_motif = ""
    for i in range(0,len(motif)):
        letter = motif[i]
        if letter in iupacdict:
            transl_motif = transl_motif + iupacdict[letter]
        else:
            # if letter not in dict use letter
            # |-- motivation: allow regEx syntax like N{3,7}
            transl_motif = transl_motif + letter
            
    return transl_motif

changeLog=\
"""Log created 2009-05-13 with the following already existing:
def fastaFileToBioSeqDict(pathToFastaFile, Alphabet='IUPACAmbiguousDNA', splitOn=None, indexAsName=0):
compl_iupacdict 
def compliment(motif, compl_iupacdict):
def reverse(text):
def revComp(seq):
def iupacList_2_regExList(motifList):
def iupac2regex(motif):

2009-05-13 -- Added def bestIdentOverLen(Bio_Blast_Record_Blast)
2009-06-17 -- Added def ifKmerInAll(kmer,dictOfSeqs, factor=1)
2009-06-23 -- Added def geneList2FastaDict(geneList, sourceFastaPath, hardMasked=True)
"""