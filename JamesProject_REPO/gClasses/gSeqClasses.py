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
        assert charsInSeq.issubset(charsInAlpha), 'Your %s sequence must contain only these characters: %s' % (self._alphabet,alphabet)

    def toString(self):    
        return ''.join(self._seq)

    def slice(self,index1,index2):
        return self.__class__(self.toString()[index1:index2])


    def calcResidueFreqs(self, customAlpha=None):
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
        return self._rFreqs[residueChar]

#!!!!
    def getResFreqTable(self):
        # will want to sort keys in _rFreqs dict and write out to a 
        # list of strings or just a formatted string
        pass


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
    _alphabet = 'DNA'
    
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

    
# !!!!!!!
    def translate(self, frame, toString=0):
        """returns PROTseq Class derived from provided frame (1,2,3,-1,-2,-3) OR a string if toString == 1
        """
        assert frame in [1,2,3,-1,-2,-3]  # frame must only be these numbers
        pass
    
    

    
    
class PROTseq(gSeq):
    """homemade class to deal with everyday Protein seq needs
    """
    import JamesDefs
    
    
    # declare that I am AA
    _alphabet = 'AA'
    
    def __init__(self,seq='', name=None, convertToUpper=1):
        gSeq.__init__(self,seq, name, convertToUpper)
        
        # enforce correct alphabet
        self._enforceAlphabet(supportVars.iupacAA)