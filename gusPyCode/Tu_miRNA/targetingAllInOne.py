from gusPyCode.defs import bioDefs
import re


class SeedModels:
    def __init__(self):
        """
        Class to centrally define the versions of seed representations and facilitate
        the separation of data pertaining to each version of the seed.

        IMPORTANT: if index[2]: should be in form:
        'instructions = ''.join([nuc_matched_in_message, pos_in_miRNA])
        """
        self.seedModels = {'m3_to_m8':[2,8],
                          'm2_to_m7':[1,7],
                          'A1_to_m7':[0,7,'A1'],
                          'm2_to_m8':[1,8],
                          'A1_to_m8':[0,8,'A1']}

class miRNA(SeedModels):
    """
    Class to represent miRNAs, including the versions of seed matches.
    """
    
    def __init__(self,miRNA_or_kmer):
        """
        Takes: mature miRNA or kmer<7bp.  If miRNA: extracts all 'versions' of seed defined in
        SeedModels and calcs Complement for matching message; stores under 'seed' and
        'seedComp' respectively.  If kmer<7bp: assumes direct dictation (does NOT derive
        multiple versions), calcs Complement for matching message; stores under 'seed'
        and 'seedComp' respectively.
        """
        
        # Inherit supported seed models.
        SeedModels.__init__(self)
        
        # Define obj vars for clairity
        self.sourceSeq     = miRNA_or_kmer.upper()
        self.sourceSeqType = None
        self.seed          = {}
        self.seedComp      = {}
        
        self._loadSeed()
        self._loadSeedComp()

    def _loadSeed(self):
        # Determine whether we got miRNA or kmer and initialize
        # 'seed' and 'seedComp' dicts.
        if len(self.sourceSeq) > 8:
            self.sourceSeqType = 'miRNA'
            # Extract all versions of seed from miRNA
            for seedType in self.seedModels:
                # If no extra info take the defined slice
                if len(self.seedModels[seedType]) == 2:
                    index_0 = self.seedModels[seedType][0]
                    index_1 = self.seedModels[seedType][1]
                    self.seed[seedType] = self.sourceSeq[index_0:index_1].replace('U','T') 
                # If extra data: use 'instructions' to construct seed where
                # revComp will generate match to message.
                elif len(self.seedModels[seedType]) == 3:
                    nuc,pos = list(self.seedModels[seedType][2])
                    pos = int(pos) 
                    index_0 = self.seedModels[seedType][0]
                    index_1 = self.seedModels[seedType][1]
                    self.seed[seedType] = self.sourceSeq[index_0:index_1].replace('U','T') 
                    # Convert correct position to reflect adjustment listed
                    # in 'instructions'.
                    self.seed[seedType] = list(self.seed[seedType])    # convert to list for inplace mutation
                    self.seed[seedType][pos-1] = bioDefs.revComp(nuc)  # insert rvCmp'd nuc
                    self.seed[seedType] = ''.join(self.seed[seedType]) # convert back to string
            
        else:
            self.sourceSeqType = 'kmer'
            self.seed['kmer'] = self.sourceSeq.replace('U','T')
    
    def _loadSeedComp(self):
        for seedType in self.seed:
            self.seedComp[seedType] = bioDefs.revComp(self.seed[seedType])

class SingleSeqHits:
    """Class to store the individual genes in each genome that contain a match to
the various versions of supported models for a seed."""
    def __init__(self, genomeTokens):
        
        # Define obj vars for clairity
        self.seedHits     = {}   # keys=str(genomeTokens) : values=dict(keys=str(geneNames) : values=Set([seedModels that succeeded]))
        self.seedCompHits = {}   # keys=str(genomeTokens) : values=dict(keys=str(geneNames) : values=Set([seedCompModels that succeeded]))
        self.genomes  = genomeTokens
        
        # Finish initializing object
        for g in self.genomes:
            self.seedHits[g]     = {}
            self.seedCompHits[g] = {}
        
class OrthologHits(SeedModels):
    """
    Class to store the orthologous genes in each genome that contain a match to
    the various versions of supported models for a seed.
    """
    def __init__(self, nOrderOrthoDefs):
        
        # Define obj vars for clairity
        self.seedHits     = {}   # keys=str(2|3|...|N) : values=dict(keys=frozenset(X001,Y033,Z584) : values=Set([seedModels that succeeded]))
        self.seedCompHits = {}   # 
        self.orthos       = nOrderOrthoDefs
        
        # Finish initializing object
        for o in self.orthos:
            self.seedHits[o]     = {}
            self.seedCompHits[o] = {}
        
 
class SeedData(SeedModels):
    """
    Class to represent the seqeunce and targeting information pertaining to a single
    seed-like sequence as well as the methods involved in collecting the information.
    """
    def __init__(self,sourceMiRNA, utrs, orthoRelations, genomeTokens):
        # Inherit supported seed models.
        SeedModels.__init__(self) 
        
        # Define obj vars for clairity
        self.genomeTokens   = genomeTokens    # List of Letter part of gene names ('AGAP')
        self.orthos         = orthoRelations  # for N genomes, dict of keys(2toN), vals(list of) 
        self.sourceMiRNA    = sourceMiRNA
        self.seedModels     = Seed(sourceSeq)   
        self.seedGeneHits   = SingleSeqHits(self.genomeTokens)  
        self.seedOrthoHits  = OrthologHits(self.orthos)
        
        
    
        
        
    #def _learnGenomes(self):
        #"""Uses self.orthos[N] to determine the genome prefixes of all genomes, then
#updates self.genomeTokens."""
        
        ## Learn largest-way combos
        #allWayKey = max(self.orthos.keys())
        ## Get list of geneSets
        #geneNames = list(self.orthos[allWayKey].keys()[0])
        ## Extract letter portion of gene tokens
        #tokens = []
        #lettersRegEx = re.compile('^\D+', re.IGNORECASE)
        #for name in geneNames:
            #assert lettersRegEx.search(name) != None , 'geneNames(%s) do not seem to have form "PREFIX000000"' % (name)
            #tokens.append(lettersRegEx.search(name).group())
            #print 'ln141: Token = %s' % (lettersRegEx.search(name).group())
        
        ## Deposit in genomeTokens
        #self.genomeTokens = tokens
        

    
    def getMatches(self, seqDict):
        pass


#-#  Global defs  #-#
def learnOrthoRelations():
    pass
def learnGenomes():
    pass
    

    


    
changeLog = """2009-05-25 -- script created.
2009-05-25 -- skeletoned some attribs/meths
2009-05-25 -- coded class SeedModels
2009-05-25 -- coded class Seed(SeedModels)
2009-06-16 -- changed Seed class name to miRNA"""

notes = """2009-05-25 -- Recording hits to each type of seed:
\t hit == dict{key=geneName : value=setOfModelsFound[model1,model2]}?
2009-05-25 -- Should move learnGenomes to MAIN so that I only do it once"""