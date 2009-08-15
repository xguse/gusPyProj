import re
import xpermutations
import JamesDefs
import bioDefs





#variable to centrally define the versions of seed representations and facilitate
#the separation of data pertaining to each version of the seed.

#IMPORTANT: if index[2]: it should be in form:
#''.join([nuc_matched_in_message, pos_in_miRNA]),   pos_in_miRNA starting with 1 NOT 0!!

_seedModels = {'m3_to_m8':[2,8],
              'm2_to_m7':[1,7],
              'A1_to_m7':[0,7,'A1'],
              'm2_to_m8':[1,8],
              'A1_to_m8':[0,8,'A1'],}

class miRNA:
    """
    Class to represent miRNAs, including the versions of seed matches.
    """
    
    def __init__(self,miRNA,set4realMatches,orthoRelations=None,genomeTokens=None,name=None,numOfCtrls=15):
        """
        Takes: mature miRNA or kmer<7bp.  If miRNA: extracts all 'versions' of seed defined in
        SeedModels and calcs corresponding mRNA matching seq and derives a series of control seqs;
        stores under 'seedMatchs' and 'seedCtrls' respectively.
        
        Ctrls are derived AFTER self.__init__ with obj method 'buildCtrls()' by randomly permuting the
        nucleotides in the seedMatch version while restricting the resulting list from containing any
        of the seedMatches recorded in provided list4realMatches during a loop of miRNA obj initiations.
        """
        

        # Define obj vars for clairity
        self.name           = name
        self.sourceSeq      = miRNA.upper()
        self.sourceSeqType  = None
        self.genomes        = genomeTokens
        self.orthos         = orthoRelations
        self.matchVersions  = {}  # formerly self.seed
        self.matchData      = {}  # keys=str(seedType) : values=set([matchHits])
        self.ctrlData       = {}  # keys=str(seedType) : values=[set([ctrlVerOneHits]),set([ctrlVerTwoHits]), ... ]
        self.matchCounts    = {}  # keys=str(seedType) : values=[countInAnyGenome,countIn1:1s, countIn1:1:1s]
        self.ctrlCounts     = {}  # keys=str(seedType) : values=[[countInAnyGenome_1,countIn1:1s_1, countIn1:1:1s_1],[countInAnyGenome_2,countIn1:1s_2, countIn1:1:1s_2] ...]
        
        self._loadmatchVersions(set4realMatches)
        
        # Data storage variables are initialized during the counting,after
        # the ctrls have been specified
        


    def _loadmatchVersions(self, set4realMatches):
        # Determine whether we got miRNA or kmer and initialize
        # matchVersions Dict.
        assert self.sourceSeq > 8, \
               'ERROR: self.sourceSeq (%s) must be > 8 nt long.' % (self.sourceSeq)
        self.sourceSeqType = 'miRNA'
        # We use the m2_to_m8 version to represent the seedVersion set for an miRNA since all versions can be constructed from this one.
        set4realMatches.add(bioDefs.revComp(self.sourceSeq[1:8].replace('U','T')))
        
        matchVersions = self._buildMatchVersions(self.sourceSeq)
        for seedType in _seedModels:
            self.matchVersions[seedType] = matchVersions[seedType]
            # Log each version of the seed match to use as restrictedList when building Ctrls
            set4realMatches.add(matchVersions[seedType][0]) # for technical reasons matchVersions[seedType] == ['aSeedMatchSeq'] so we must index it

    # =========
    
    def _buildMatchVersions(self, miRNA_seq):
        """
        
        """
        returnDict = {}
        # Extract all versions of seed from miRNA
        for seedType in _seedModels:
            returnDict[seedType] = []
            # If NO third index like here: [0,7,'A1']
            if len(_seedModels[seedType]) == 2:
                index_0 = _seedModels[seedType][0]
                index_1 = _seedModels[seedType][1]
                returnDict[seedType].append(bioDefs.revComp(miRNA_seq[index_0:index_1].replace('U','T'))) 
            # If extra data: use 'instructions index' to construct seedVersion
            elif len(_seedModels[seedType]) == 3:
                nuc,pos = list(_seedModels[seedType][2])
                pos = int(pos) 
                index_0 = _seedModels[seedType][0]
                index_1 = _seedModels[seedType][1]
                returnDict[seedType].append(bioDefs.revComp(miRNA_seq[index_0:index_1].replace('U','T'))) 
                # Convert correct position to reflect adjustment listed
                # in 'instructions'.
                returnDict[seedType][0]       = list(returnDict[seedType][0])    # convert to list for inplace mutation
                returnDict[seedType][0][-pos] = nuc                              # insert nuc
                returnDict[seedType][0]       = ''.join(returnDict[seedType][0]) # convert back to string
        
        return returnDict
    
    def buildCtrlsFromMatchVers(self, restrictedList, numOfCtrls=15):
        """
        WARNING: This should only be called after the entire list of real seeds have been initialized!
        
        Computes 15 permutations of each 'true' matchVersion and screens them for real seqs in restrictedList
        to prevent using KNOWN seed matches for controls. Ctrl seed matches are stored in a list located
        at second index of the list located in dict entry self.matchVersions[seedType]. 
        """
        
        # check to see whether this has already been done.
        # If so, complain and die.
        # Else, append an empty list as index_1 after REAL matchSeq for each version
        for seedType in self.matchVersions:
            assert type(self.matchVersions[seedType]) == type([]), \
                   'ERROR: %s.matchVersions[%s] is not type: list.' % (self.name,seedType)
            assert len(self.matchVersions[seedType]) == 1, \
                   'ERROR: len(%s.matchVersions[%s]) is not 1; ctrls seqs may have already been built.' % (self.name,seedType)
            self.matchVersions[seedType].append([])
        
        # permute and screen each seedVersion
        for seedType in self.matchVersions:
            # If NO third index like here: [0,7,'A1']
            if len(_seedModels[seedType]) == 2:
                # Select 15 random permutations of matchVersions[seedType] that are not in the 
                # restrictedList.
                matchPermList = [''.join(x) for x in xpermutations.xpermutations(list(self.matchVersions[seedType][0]))]
                while len(self.matchVersions[seedType][1]) < numOfCtrls:
                    permSeq = JamesDefs.randFromList_noReplace(matchPermList)
                    if permSeq not in restrictedList:
                        # Append permuted Seq if not in restrictedList
                        self.matchVersions[seedType][1].append(permSeq)
            # If extra data: use 'instructions index' to only permute the nucs not explicitly
            # defined in the seedModel
            elif len(_seedModels[seedType]) == 3:
                nuc,pos = list(_seedModels[seedType][2])
                # Leave 1-registered bc we will use negIndex bc dealing with rvCmp of miRNA
                # so pos == 1 actually means pos == LAST in matchSeq 
                pos = int(pos) 
                # explode seq to remove defined nuc in place
                seq2Perm = list(self.matchVersions[seedType][0])
                del seq2Perm[-pos]
                # Generate permutations from remaining nucs,
                matchPermList = [x for x in xpermutations.xpermutations(seq2Perm)]
                while len(self.matchVersions[seedType][1]) < numOfCtrls:
                    permSeq = JamesDefs.randFromList_noReplace(matchPermList)
                    # Replace nuc and check restricted list.
                    if pos > 1: permSeq.insert(-pos+1,nuc)
                    else:       permSeq.append(nuc)
                    permSeq = ''.join(permSeq)
                    if permSeq not in restrictedList:
                        # Append permuted Seq if not in restrictedList
                        self.matchVersions[seedType][1].append(permSeq) 
            

    def buildCtrlsFromProSeed(self, restrictedList, numOfCtrls=15):
        """
        WARNING: This should only be called after the entire list of real seeds have been initialized!
        
        Computes 15 permutations of the 'true' proSeed matching seqeunce (m2_to_m8) and derives
        matchVersions as in the true case. The permuted sequence is checked agaisnt the restrictedList
        to prevent using KNOWN seed matches for controls. Ctrl seed matches are stored in a list located
        at second index of the list located in dict entry self.matchVersions[seedType]. Each seedVersion 
        of a ctrl set will share the same index number.
        """
        
        # check to see whether this has already been done.
        # If so, complain and die.
        # Else, append an empty list as index_1 after REAL matchSeq for each version
        for seedType in self.matchVersions:
            assert type(self.matchVersions[seedType]) == type([]), \
                   'ERROR: %s.matchVersions[%s] is not type: list.' % (self.name,seedType)
            assert len(self.matchVersions[seedType]) == 1, \
                   'ERROR: len(%s.matchVersions[%s]) is not 1; ctrls seqs may have already been built.' % (self.name,seedType)
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
        
        # Use each chosenSeq to generate the diff matchVersions
        for seq in chosenPerms:
            # Create Fake miRNA with seq at the seed location to feed to _buildMatchVersions()
            seq = 'N%sNNNNNNNNNNNNN' % (bioDefs.revComp(seq))
            matchVersions = self._buildMatchVersions(seq)
            for seedType in self.matchVersions:
                self.matchVersions[seedType][1].append(matchVersions[seedType][0]) # must use index[0] bc _buildMatchVersions returns a list len==1 
            
    def tallyHits(self, flankingRegions):
        """
        Takes a dict of flanking DNA seqs.  For each seq, it updates the counts for each real
        matchVersion and each permutation of the ctrls for each matchVersion in self.matchData,
        and self.ctrlData.
        """
## ! ## ! ## ! ## !  ##  MUST STILL RUN TEST FOR THIS FUNCTION  ## ! ## ! ## ! ## ! ##
        
        # Initialize data variables
        for seedType in _seedModels:
            self.matchData[seedType] = set()
            self.ctrlData[seedType]  = [set()]*len(self.matchVersions[seedType][1])
        
        # Cycle through flanking regions
        for gene in flankingRegions:
            # Cycle through seedTypes
            for seedType in _seedModels:
                # Record hits for matchVersion
                if flankingRegions[gene].find(self.matchVersions[seedType][0]) >= 0:
                    self.matchData[seedType].add(gene)
                # Record hits for each ctrl permutation of the matchVersion
                for i in range(len(self.matchVersions[seedType][1])):
                    if flankingRegions[gene].find(self.matchVersions[seedType][1][i]) >= 0:
                        self.ctrlData[seedType][i].add(gene)


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
2009-08-07 -- copied to this file from targetingAllInOne
2009-08-13 -- added miRNA.buildCtrlsFromMatchVers() and renamed miRNA.buildCtrls() to buildCtrlsFromProSeed()
2009-08-13 -- changed miRNA.__init__() to log all matchVersions to use for the restrictedList
2009-08-13 --  altered var names of class SingleSeqHits
2009-08-14 -- declassified SeedModels in favor of a global variable '_seedModels' (the class had no methods and I did not plan for it to have any)
2009-08-15 -- Did away with classes SingleHits and OrthoHits, moved data storage into miRNA
2009-08-15 -- Wrote miRNA.tallyHits
"""

notes = """2009-05-25 -- Recording hits to each type of seed:
\t hit == dict{key=geneName : value=setOfModelsFound[model1,model2]}?
2009-05-25 -- Should move learnGenomes to MAIN so that I only do it once"""