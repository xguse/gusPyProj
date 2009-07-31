


class OrthoGroup:
    """Sets up an orthoGroup object when given a Dict of form: KEYS:<geneNames>, VALUES:<DNAsequence>
    Manages the searching of self with TAMO motif object, and modifies 'hitList' attribute in self
    to reflect k-mer searched and positive fraction."""
    


    
    
    
    def __init__(self, orthoDict):
        
        # object data 
        self.hitList = []
        self.orthoDict = orthoDict
        
        self._validate()
        
    def _validate(self):
        assert type(self.orthoDict) == type({}), \
               '''orthoSet.orthoDict must be a dictionary.
               You provided:%s which was %s'''\
               % (self.orthoDict, type(self.orthoDict))
        
    
    def searchMotif(self,TAMO_motif_obj,scoreFactor=0.75):
        """Scans all seqs in self using TAMO motif obj and scoreFactor.  Records k-mer string and fraction
        of positive seqs in self into self.hitlist."""
        
        assert type(scoreFactor) == type(1) or type(scoreFactor) == type(0.1), \
               '''orthoSet.searchMotif argument "scoreFactor" must be of type int or float.
               You provided: %s which was %s'''\
               % (scoreFactor,type(scoreFactor))
        
        hits = 0
        for seq in self.orthoDict:
            if TAMO_motif_obj.scan(self.orthoDict[seq],factor=scoreFactor)[0]: # ALWAYS returns a tuple, so check first index
                hits+=1
        
        # update hitList with: [kMerString, fracOfPositives, scoreFactorUsed]
        self.hitList.append([TAMO_motif_obj.oneletter, hits/float(len(self.orthoDict)), scoreFactor])






changeLog = \
"""2009-03-28 -- creation
2009-03-28 -- added orthoSet.__init__()
2009-03-28 -- added orthoSet.searchMotif()
2009-03-28 -- Changed orthoSet class name to 'OrthoGroup' bc it is not a PySetObj and I tend to use that convention

"""


