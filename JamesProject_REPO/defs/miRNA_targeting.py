import re
import xpermutations
import JamesDefs
import bioDefs



class SeedModels:
    def __init__(self):
        """
        Class to centrally define the versions of seed representations and facilitate
        the separation of data pertaining to each version of the seed.

        IMPORTANT: if index[2]: it should be in form:
        ''.join([nuc_matched_in_message, pos_in_miRNA]),   pos_in_miRNA starting with 1 NOT 0!!
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
    
    def __init__(self,miRNA_or_kmer,set4realMatches, name=None):
        # will probably remove the option for submitting a kmer
        """
        Takes: mature miRNA or kmer<7bp.  If miRNA: extracts all 'versions' of seed defined in
        SeedModels and calcs corresponding mRNA matching seq and derives a series of control seqs;
        stores under 'seedMatchs' and 'seedCtrls' respectively.  If kmer<7bp: assumes direct dictation
        (does NOT derive multiple versions), corresponding mRNA matching seq and derives a series of
        control seqs; stores under 'seedMatchs' and 'seedCtrls' respectively.
        
        Ctrls are derived AFTER self.__init__ with obj method 'buildCtrls()' by randomly permuting the
        nucleotides in the seedMatch version while restricting the resulting list from containing any
        of the seedMatches recorded in provided list4realMatches during a loop of miRNA obj initiations.
        """
        
        # Inherit supported seed models.
        SeedModels.__init__(self)
        
        # Define obj vars for clairity
        self.name          = name
        self.sourceSeq     = miRNA_or_kmer.upper()
        self.sourceSeqType = None
        self.matchVersions  = {}  # formerly self.seed
        
        self._loadmatchVersions(set4realMatches)


    def _loadmatchVersions(self, set4realMatches):
        # Determine whether we got miRNA or kmer and initialize
        # matchVersions Dict.
        assert self.sourceSeq > 8, \
               'ERROR: self.sourceSeq (%s) must be > 8 nt long.' % (self.sourceSeq)
        self.sourceSeqType = 'miRNA'
        # We use the m2_to_m8 version to represent the seedVersion set for an miRNA since all versions can be constructed from this one.
        set4realMatches.add(bioDefs.revComp(self.sourceSeq[1:8].replace('U','T')))
        
        matchVersions = self._buildMatchVersions(self.sourceSeq)
        for seedType in self.seedModels:
            self.matchVersions[seedType] = matchVersions[seedType]

    # =========
    
    def _buildMatchVersions(self, miRNA_seq):
        """
        
        """
        returnDict = {}
        # Extract all versions of seed from miRNA
        for seedType in self.seedModels:
            returnDict[seedType] = []
            # If NO third index like here: [0,7,'A1']
            if len(self.seedModels[seedType]) == 2:
                index_0 = self.seedModels[seedType][0]
                index_1 = self.seedModels[seedType][1]
                returnDict[seedType].append(bioDefs.revComp(miRNA_seq[index_0:index_1].replace('U','T'))) 
            # If extra data: use 'instructions index' to construct seedVersion
            elif len(self.seedModels[seedType]) == 3:
                nuc,pos = list(self.seedModels[seedType][2])
                pos = int(pos) 
                index_0 = self.seedModels[seedType][0]
                index_1 = self.seedModels[seedType][1]
                returnDict[seedType].append(bioDefs.revComp(miRNA_seq[index_0:index_1].replace('U','T'))) 
                # Convert correct position to reflect adjustment listed
                # in 'instructions'.
                returnDict[seedType][0]       = list(returnDict[seedType][0])    # convert to list for inplace mutation
                returnDict[seedType][0][-pos] = nuc                              # insert nuc
                returnDict[seedType][0]       = ''.join(returnDict[seedType][0]) # convert back to string
        
        return returnDict
    
    def buildCtrls(self, restrictedList, numOfCtrls=15):
        """
        WARNING: This should only be called after the entire list of real seeds have been initialized!
        
        Computes 15 permutations of the 'true' seed matching seqeunce and derives matchVersions
        as in the true case. The permuted sequence is checked agaisnt the resttictedList to prevent
        using KNOWN seed matches for controls. Ctrl seed matches are stored in a list located at
        second index of the list located in dict entry self.matchVersions[seedType]. Each seedVersion 
        of a ctrl set will share the same index number.
        """
        
        # check to see whether this has already been done.
        # If so, complain and die.
        # Else, append an empty list as index_1 after REAL matchSeq for each version
        for seedType in self.matchVersions:
            assert type(self.matchVersions[seedType]) == type([]), \
                   'ERROR: %s.matchVersions[%s] is not type: list.' % (self.name,seedType)
            assert len(self.matchVersions[seedType]) == 1, \
                   'ERROR: len(%s.matchVersions[%s]) is not 1; buildCtrls() may have already been run.' % (self.name,seedType)
            self.matchVersions[seedType].append([])
        
        proSeed = self.sourceSeq[1:8]
        matchPerms  = [''.join(x) for x in xpermutations.xpermutations(list(self.matchVersions['m2_to_m8'][0]))]
        
        # Select 15 random permutations of matchVersions['m2_to_m8'] that are not in the 
        # restrictedList.
        chosenPerms = []
        while len(chosenPerms) < numOfCtrls:
            permSeq = JamesDefs.randFromList_noReplace(matchPerms)
            if permSeq not in restrictedList:
                chosenPerms.append(permSeq)
        
        # this list will be ~40K long.  Lets explicitly kill it for memory's sake
        del matchPerms
        # Use each chosenSeq to generate the diff matchVersions
        for seq in chosenPerms:
            # Create Fake miRNA with seq at the seed location to feed to _buildMatchVersions()
            seq = 'N%sNNNNNNNNNNNNN' % (bioDefs.revComp(seq))
            matchVersions = self._buildMatchVersions(seq)
            for seedType in self.matchVersions:
                self.matchVersions[seedType][1].append(matchVersions[seedType][0]) # must use index[0] bc _buildMatchVersions returns a list len==1 
            

class SingleSeqHits:
    """
    Class to store the individual genes in each genome that contain a match to
    the various versions of supported models for a seed.
    """
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
2009-06-16 -- changed Seed class name to miRNA
2009-09-07 -- copied to this file from targetingAllInOne"""

notes = """2009-05-25 -- Recording hits to each type of seed:
\t hit == dict{key=geneName : value=setOfModelsFound[model1,model2]}?
2009-05-25 -- Should move learnGenomes to MAIN so that I only do it once"""