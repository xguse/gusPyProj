from __future__ import division
import sys
import supportVars 
import JamesDefs
import sets
import re
from decimal import Decimal


##class list1(list):
    ##def __init__(self,l=[],name=None):
        ##assert (type(l) == type([]))
        ##list.__init__(self,l)
        ##self.name = name

class gSeq:
    """Parent class for Seq objects
    """
    
    
    # Default Variables
    _alphabet = None
    # initiate a Dict to store int's rFrequencies
    _rFreqs = {}
    
    

    def __init__(self,seq='', name=None, convertToUpper=1):        
        # enforce string type data
        assert (type(seq) == type("") or type(seq) == type(u""))  # must use a string but can be a unicode string
        
        if convertToUpper != 0:
            self._seq = list(seq.upper())
        else:
            self._seq = list(seq)
        self.name = name
        
        
        

    def __add__(self, otherSeqObj):
        """returns new seqObj with seq2 appened to end of seq1
        """
        # enforce same class type
        assert self.__class__.__name__ == otherSeqObj.__class__.__name__ # can not add diff seqTypes!!
        
        # spawn new seqObj after concatentating strings
        newSeqObj = self.__class__(self.toString()+otherSeqObj.toString())
        return newSeqObj
    
    def __len__(self):
        return len(self._seq)
    
    def __repr__(self):
        return "%s(%s, %s)" % (self.__class__.__name__,
                               repr(self.toString()),
                               repr(self._alphabet))
    def __str__(self):
        if len(self._seq) > 60:
            s = repr(self._seq[:60] + " ...")
        else:
            s = repr(self._seq)
        return "%s(%s, %s)" % (self.__class__.__name__, s,
                               repr(self._alphabet))
    
    def _enforceAlphabet(self,alphabet):
        # enforce correct alphabet
        charsInSeq   = sets.Set(self._seq)
        charsInAlpha = sets.Set(alphabet)

#!!     will want to build in a way to report the offending Char (prob use Set methods)
        assert charsInSeq.issubset(charsInAlpha), 'Your %s sequence contains %s. The following are allowed: %s' % (self._alphabet,list(charsInSeq.difference(charsInAlpha)),alphabet)

    def toString(self):    
        return ''.join(self._seq)

    def slice(self,index1,index2):
        return self.__class__(self.toString()[index1:index2])


    def _calcResidueFreqs(self, customAlpha=None):
        """Usage: residueFreqs(customAlpha=None) customAlpha is list of characters
        Returns a Dict with Key<residue>,Value<Freq of Key>"""
        # set correct alphabet to use
        if type(customAlpha) is list:
            alpha = customAlpha
        elif self._alphabet == "DNA":
            alpha = supportVars.iupacDNA
        elif self._alphabet == "AA":
            alpha = supportVars.iupacAA
        
        for char in alpha:
            self._rFreqs[char] = Decimal(self._seq.count(char))/len(self._seq)
    
    def getResFreq(self, residueChar):
        assert residueChar in supportVars.self._alphabet
        ##if self._rFreqs:
            ##return self._rFreqs[residueChar]
        ##else:
        self._calcResidueFreqs()
        return self._rFreqs[residueChar]

#!!!!
    def getResFreqTable(self):
        # will want to sort keys in _rFreqs dict and write out to a 
        # list of strings or just a formatted string
        self._calcResidueFreqs()
        sortedResidues = self._rFreqs.keys()
        sortedResidues.sort()
        
        freqTab = ['Residue\t%s' % (self.name)]
        for res in sortedResidues:
            freqTab.append(self._rFreqs[res])
        return freqTab


# !!!!!    
    def search(self,):
        assert regEx == 0 or 1  # 0=unambiguous pattern match; 1=use regular expressions
        if regEx == 1:
            pass
        else:
            pass
        
    def toFasta(self):
        return ">%s\n%s\n" % (self.name, self.toString())
        
        


class DNAseq(gSeq):
    """homemade class to deal with everyday DNA seq needs
    """
    import JamesDefs
    
    
    # declare that I am DNA
    _alphabet = 'iupacDNA'
    
    def __init__(self,seq='', name=None, convertToUpper=1):
        gSeq.__init__(self,seq, name, convertToUpper)
        
        # enforce correct alphabet
        self._enforceAlphabet(supportVars.iupacDNA)



    def revCmp(self, toString=0):
        assert toString == 0 or 1 # 0 = return new SeqObj; 1 = return string
        
        if toString == 1:
            return JamesDefs.revComp(self.toString())
        else:
            return self.__class__(JamesDefs.revComp(self.toString()))

    def codons(self, frame=1): 
        """Return list of codons for self._seq""" 
        assert frame in [1,2,3,-1,-2,-3] , 'Frame must only be these numbers.'
        if frame < 0:
            rvCmp = self.revCmp()
            end = len(rvCmp._seq) - (len(rvCmp._seq) % 3) - 1 
            codons = [''.join(rvCmp._seq[i:i+3]) for i in range(0+abs(frame)-1, end, 3)] 
            return codons 
        else:
            end = len(self._seq) - (len(self._seq) % 3) - 1 
            codons = [''.join(self._seq[i:i+3]) for i in range(0+frame-1, end, 3)] 
            return codons 
    
# !!!!!!!
    def translate(self, frame=1, toString=0):
        """returns PROTseq Class derived from provided frame (1,2,3,-1,-2,-3) OR a string if toString == 1
        """
        assert frame in [1,2,3,-1,-2,-3] , 'Frame must only be these numbers.'
        
        # get list of codons
        codons = self.codons(frame)
        
        peptide = ''
        for each in codons:
            if len(each) == 3:
                peptide =  peptide+supportVars.geneCode[each]
                
        return PROTseq(peptide)
    
    

    
    
class PROTseq(gSeq):
    """homemade class to deal with everyday Protein seq needs
    """
    import JamesDefs
    
    
    # declare that I am AA
    _alphabet = 'iupacAA'
    
    def __init__(self,seq='', name=None, convertToUpper=1):
        gSeq.__init__(self,seq, name, convertToUpper)
        
        # enforce correct alphabet
        self._enforceAlphabet(supportVars.iupacAA)
        
        
        
        
        
        
################ Change Log ######################
##gSeq:

##DNAseq:
##- added codons()
##- completed version 1 of translate() 
##PROTseq: