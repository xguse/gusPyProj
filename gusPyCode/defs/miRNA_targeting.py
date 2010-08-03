import numpy
import re
from random import sample
from TAMO.seq import Fasta
from gusPyCode.defs import xpermutations
from gusPyCode.defs import JamesDefs
from gusPyCode.defs import bioDefs
from gusPyCode.defs import mathDefs
from pprint import pprint


def shuffSeedsInMatMiRs(inPath,outPathBase,seed=[1,8],num=1):
    """Takes mature miR fasta file from inPath, permutes nucs in <seed>
    <num> times and writes records to outPaths derived from outPathBase
    with unique ctrl numbers.
    
    <seed> coords are as pyhton slice coords (0-based and
    non-inclusive on the right side)"""
    inFile = open(inPath,'rU')
    outFiles = []
    for i in range(num):
        outFiles.append(open('%s.ctrl%s.fas' % (outPathBase,i),'w'))

    for line in inFile:
        if line.startswith('>'):
            for oF in outFiles:
                oF.write(line)
        else:
            line = line.strip('\n')
            seedSeq  = line[seed[0]:seed[1]]
            fvPrm    = line[:seed[0]]
            trPrm    = line[seed[1]:]
            shfSeeds = [x for x in xpermutations.xpermutations(list(seedSeq))]
            sampShfs = sample(shfSeeds,num)
            for i in range(num):
                nLine = '%s%s%s\n' % (fvPrm,''.join(sampShfs[i]),trPrm)
                outFiles[i].write(nLine)
            


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
    Class to represent seedMatches.
    """
    
    def __init__(self,miRNA,set4realMatches,orthoRelations=None,name=None,numOfCtrls=15):
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
        ##self.genomes        = genomeTokens
        self.orthos         = orthoRelations
        self.matchVersions  = {}  # formerly self.seed
        self.matchData      = {}  # keys=str(seedType) : values=set([matchHits])
        self.ctrlData       = {}  # keys=str(seedType) : values=[set([ctrlVerOneHits]),set([ctrlVerTwoHits]), ... ]
        self.matchEvents    = {}  # keys=str(seedType) : values=[None,genesInAnyGenome,pairsInAtLeast1:1s, tripsIn1:1:1s]
        self.ctrlEvents     = {}  # keys=str(seedType) : values=[[None,genesInAnyGenome_1,pairsInAtLeast1:1s_1, tripsIn1:1:1s_1],[None,genesInAnyGenome_2,pairsInAtLeast1:1s_2, tripsIn1:1:1s_2]...]
        self.matchCounts    = {}  # keys=str(seedType) : values=[countInNone,countInAnyGenome,countInAtLeast1:1s, countIn1:1:1s]
        self.ctrlCounts     = {}  # keys=str(seedType) : values=[[countInNone_1, countInAnyGenome_1,countInAtLeast1:1s_1, countIn1:1:1s_1],[countInNone_2, countInAnyGenome_2,countInAtLeast1:1s_2, countIn1:1:1s_2] ...]
        self.ctrlMeanStd    = {}  # keys=str(seedType) : values=[[countInNone-mean, countInNone-stdv], [countInAnyGenome-mean, countInAnyGenome-stdv],[countInAtLeast1:1s-mean, countInAtLeast1:1s-stdv], [countIn1:1:1s-mean, countIn1:1:1s-stdv]]
        self.ctrlMedianStd  = {}  # keys=str(seedType) : values=[[countInNone-median, countInNone-stdv], [countInAnyGenome-median, countInAnyGenome-stdv],[countInAtLeast1:1s-median, countInAtLeast1:1s-stdv], [countIn1:1:1s-median, countIn1:1:1s-stdv]]
        self.ctrlMedianMAD  = {}  # keys=str(seedType) : values=[[countInNone-median, countInNone-MAD], [countInAnyGenome-median, countInAnyGenome-MAD],[countInAtLeast1:1s-median, countInAtLeast1:1s-MAD], [countIn1:1:1s-median, countIn1:1:1s-MAD]]
        self.seedScores     = {}  # keys=str(seedType) : values={keys=scoreType:values=zScoreNumber}
        
        self._loadmatchVersions(set4realMatches)
        
        # Structures of data storage variables are initialized in the def that fills them in.
        


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
        
        Computes usr supplied number of permutations (dflt=15) of each 'true' matchVersion and screens them
        for real seqs in restrictedList to prevent using KNOWN seed matches for controls. Ctrl seed matches
        are stored in a list located at second index of the list located in dict entry
        self.matchVersions[seedType]. 
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
                # Select numOfCtrls random permutations of matchVersions[seedType] that are not in the 
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
        ##assert True==False, \
               ##"""WARNING!!! miRNA.buildCtrlsFromMatchVers() should be used instead!!
               ##If you REALLY want to use this method, modify or remove this assert statement.
               
               ##But seriously...  use miRNA.buildCtrlsFromMatchVers()."""
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
        Takes a dict of flanking DNA seqs.  For each seq, it updates the genes containing each
        real matchVersion and each permutation of the ctrls for each matchVersion in self.matchData,
        and self.ctrlData.
        """
        
        # Initialize data variables
        for seedType in _seedModels:
            self.matchData[seedType] = set()
            self.ctrlData[seedType]  = []
            for i in range(len(self.matchVersions[seedType][1])):
                self.ctrlData[seedType].append(set())
        
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

    def countHitsInOrthos(self,returnGenes=True):
        """
        Uses results of miRNA.tallyHits() and self.orthos to count how many genes the miRNA seed
        hits in at least one genome, in at least two orthologs, and in all three orthologs.  If
        returnGenes: returns tuple of two dicts:
        matchDict(keys=seedType : vals=[None,genesWithMatch,genePairsWithMatch,geneTriplesWithMatch])
        ctrlDict(keys=seedType : vals=[[],genesWithMatch_1,genePairsWithMatch_1,geneTriplesWithMatch_1],
                                       [],genesWithMatch_2,genePairsWithMatch_2,geneTriplesWithMatch_2],
                                       ...])
        """
        # make sure we have tallied the hits already.
        assert self.matchData and self.ctrlData, \
               'ERROR:  It looks like we have not tallied the hits yet. Call miRNA.tallyHits() first.'
        
        if returnGenes:
            rGeneNames = {}
            rCtrlNames = {}
            for seedType in _seedModels:
                rGeneNames[seedType] = [None,[],[],[]]
                rCtrlNames[seedType] = JamesDefs.initList(len(self.matchVersions[seedType][1]),[None,[],[],[]])
             
        # Initialize self.matchCounts/self.ctrlCounts
        for seedType in _seedModels:
            self.matchCounts[seedType] = [0,0,0,0]
            self.ctrlCounts[seedType]  = [[0]*4 for i in range(len(self.matchVersions[seedType][1]))]
        # Cycle through self.orthos 
        for orthoSet in self.orthos:
            assert len(orthoSet) == 3,\
                   'ERROR: It seems len(%s) != 3.'
            # Query the matcheData and ctrlData for hits in orthoSet
            for seedType in _seedModels:
                genesInMatchD = 0
                genesInCtrlD  = [0]*len(self.matchVersions[seedType][1])
                if returnGenes:
                    geneNames = [] 
                    ctrlNames = JamesDefs.initList(len(self.matchVersions[seedType][1]),[])
                # Count how many genes in each orthoSet were hit by the respective seedTypes
                for gene in orthoSet:
                    if gene in self.matchData[seedType]:
                        genesInMatchD += 1
                        if returnGenes: geneNames.append(gene)
                    for i in range(len(self.ctrlData[seedType])):
                        if gene in self.ctrlData[seedType][i]:
                            genesInCtrlD[i] += 1
                            if returnGenes: ctrlNames[i].append(gene)
                        
                # Update self.matchData based on how many hits the orthoSet got for seedType
                if genesInMatchD == 0:
                    self.matchCounts[seedType][0] += 3
                elif genesInMatchD == 1:
                    self.matchCounts[seedType][0] += 2
                    self.matchCounts[seedType][1] += 1
                    if returnGenes:
                        rGeneNames[seedType][1].extend(geneNames)
                elif genesInMatchD == 2:
                    self.matchCounts[seedType][0] += 1
                    self.matchCounts[seedType][1] += 2
                    self.matchCounts[seedType][2] += 1
                    if returnGenes:
                        rGeneNames[seedType][1].extend(geneNames)
                        rGeneNames[seedType][2].append(tuple(sorted(geneNames)))
                elif genesInMatchD == 3:
                    self.matchCounts[seedType][1] += 3
                    self.matchCounts[seedType][2] += 3
                    self.matchCounts[seedType][3] += 1
                    if returnGenes:
                        rGeneNames[seedType][1].extend(geneNames)
                        type2 = [tuple(sorted(x)) for x in xpermutations.xuniqueCombinations(geneNames,2)]
                        type3 = tuple(sorted(geneNames))
                        rGeneNames[seedType][2].extend(type2)
                        rGeneNames[seedType][3].append(type3)
                # Update self.ctrlData based on how many hits the orthoSet got in each ctrl for seedType
                for i in range(len(self.ctrlData[seedType])):
                    if genesInCtrlD[i] == 0:
                        self.ctrlCounts[seedType][i][0] += 3
                    elif genesInCtrlD[i] == 1: 
                        self.ctrlCounts[seedType][i][1] += 1
                        if returnGenes:
                            rCtrlNames[seedType][i][1].extend(ctrlNames[i])
                    elif genesInCtrlD[i] == 2:
                        self.ctrlCounts[seedType][i][1] += 2
                        self.ctrlCounts[seedType][i][2] += 1
                        if returnGenes:
                            rCtrlNames[seedType][i][1].extend(ctrlNames[i])
                            rCtrlNames[seedType][i][2].append(tuple(sorted(ctrlNames[i])))
                    elif genesInCtrlD[i] == 3:
                        self.ctrlCounts[seedType][i][1] += 3
                        self.ctrlCounts[seedType][i][2] += 3
                        self.ctrlCounts[seedType][i][3] += 1
                        if returnGenes:
                            rCtrlNames[seedType][i][1].extend(ctrlNames[i])
                            type2 = [tuple(sorted(x)) for x in xpermutations.xuniqueCombinations(ctrlNames[i],2)]
                            type3 = tuple(sorted(ctrlNames[i]))
                            rCtrlNames[seedType][i][2].extend(type2)
                            rCtrlNames[seedType][i][3].append(type3)

                if returnGenes:
                    for i in range(1,4):
                        assert len(rGeneNames[seedType][i]) == len(set(rGeneNames[seedType][i])),\
                               "ERROR: rGeneNames[%s] in miRNA(%s) has redundancy." % (i, self.name)
        if returnGenes:
            # store and return rGeneNames
            self.matchEvents = rGeneNames
            self.ctrlEvents  = rCtrlNames
            return (rGeneNames,rCtrlNames)
    
    def countHitsInOrthos4(self,genomeToken,returnGenes=True):
        """
        Uses results of miRNA.tallyHits() and self.orthos to count how many genes the miRNA seed
        hits in at least one genome, in at least two orthologs, and in all three orthologs.  If
        returnGenes: returns tuple of two dicts:
        matchDict(keys=seedType : vals=[None,genesWithMatch,genePairsWithMatch,geneTriplesWithMatch])
        ctrlDict(keys=seedType : vals=[[],genesWithMatch_1,genePairsWithMatch_1,geneTriplesWithMatch_1],
                                       [],genesWithMatch_2,genePairsWithMatch_2,geneTriplesWithMatch_2],
                                       ...])
        """
        # make sure we have tallied the hits already.
        assert self.matchData and self.ctrlData, \
               'ERROR:  It looks like we have not tallied the hits yet. Call miRNA.tallyHits() first.'
        
        if returnGenes:
            rGeneNames = {}
            rCtrlNames = {}
            for seedType in _seedModels:
                rGeneNames[seedType] = [None,[],[],[]]
                rCtrlNames[seedType] = JamesDefs.initList(len(self.matchVersions[seedType][1]),[None,[],[],[]])
             
        # Initialize self.matchCounts/self.ctrlCounts
        for seedType in _seedModels:
            self.matchCounts[seedType] = [0,0,0,0]
            self.ctrlCounts[seedType]  = [[0]*4 for i in range(len(self.matchVersions[seedType][1]))]
        # Cycle through self.orthos 
        for orthoSet in self.orthos:
            assert len(orthoSet) == 3,\
                   'ERROR: It seems len(%s) != 3.'
            # Query the matcheData and ctrlData for hits in orthoSet
            for seedType in _seedModels:
                genesInMatchD = 0
                genesInCtrlD  = [0]*len(self.matchVersions[seedType][1])
                if returnGenes:
                    geneNames = [] 
                    ctrlNames = JamesDefs.initList(len(self.matchVersions[seedType][1]),[])
                # Count how many genes in each orthoSet were hit by the respective seedTypes
                for gene in orthoSet:
                    if gene in self.matchData[seedType]:
                        genesInMatchD += 1
                        if returnGenes: geneNames.append(gene)
                    for i in range(len(self.ctrlData[seedType])):
                        if gene in self.ctrlData[seedType][i]:
                            genesInCtrlD[i] += 1
                            if returnGenes: ctrlNames[i].append(gene)
                        
                # Update self.matchData based on how many hits the orthoSet got for seedType
                if genesInMatchD == 0:
                    ##self.matchCounts[seedType][0] += 3
                    pass
                elif genesInMatchD == 1:
                    if ''.join(geneNames).find(genomeToken) != -1:
                        ##self.matchCounts[seedType][0] += 2
                        self.matchCounts[seedType][1] += 1
                        if returnGenes:
                            rGeneNames[seedType][1].extend(geneNames)
                elif genesInMatchD == 2:
                    if ''.join(geneNames).find(genomeToken) != -1:
                        ##self.matchCounts[seedType][0] += 1
                        ##self.matchCounts[seedType][1] += 1 # self.matchCounts[seedType][1] += 2
                        ##self.matchCounts[seedType][2] += 1
                        if returnGenes:
                            rGeneNames[seedType][1].extend([x for x in geneNames if x.find(genomeToken) != -1])
                            rGeneNames[seedType][2].append(tuple(sorted(geneNames)))
                elif genesInMatchD == 3:
                    if ''.join(geneNames).find(genomeToken) != -1:
                        ##self.matchCounts[seedType][1] += 1 # self.matchCounts[seedType][1] += 3
                        ##self.matchCounts[seedType][2] += 2 # self.matchCounts[seedType][2] += 3
                        ##self.matchCounts[seedType][3] += 1
                        if returnGenes:
                            rGeneNames[seedType][1].extend([x for x in geneNames if x.find(genomeToken) != -1])
                            type2 = [tuple(sorted(x)) for x in xpermutations.xuniqueCombinations(geneNames,2) if ''.join(x).find(genomeToken) != -1]
                            type3 = tuple(sorted(geneNames))
                            rGeneNames[seedType][2].extend(type2)
                            rGeneNames[seedType][3].append(type3)
                # Update self.ctrlData based on how many hits the orthoSet got in each ctrl for seedType
                for i in range(len(self.ctrlData[seedType])):
                    if genesInCtrlD[i] == 0:
                        ##self.ctrlCounts[seedType][i][0] += 3
                        pass
                    elif genesInCtrlD[i] == 1:
                        if ''.join(ctrlNames[i]).find(genomeToken) != -1:
                            ##self.ctrlCounts[seedType][i][1] += 1
                            if returnGenes:
                                rCtrlNames[seedType][i][1].extend(ctrlNames[i])
                    elif genesInCtrlD[i] == 2:
                        if ''.join(ctrlNames[i]).find(genomeToken) != -1:
                            ##self.ctrlCounts[seedType][i][1] += 1 # self.ctrlCounts[seedType][i][1] += 2
                            ##self.ctrlCounts[seedType][i][2] += 1
                            if returnGenes:
                                rCtrlNames[seedType][i][1].extend([x for x in ctrlNames[i] if x.find(genomeToken) != -1])
                                rCtrlNames[seedType][i][2].append(tuple(sorted(ctrlNames[i])))
                    elif genesInCtrlD[i] == 3:
                        if ''.join(ctrlNames[i]).find(genomeToken) != -1:
                            ##self.ctrlCounts[seedType][i][1] += 1 # self.ctrlCounts[seedType][i][1] += 3
                            ##self.ctrlCounts[seedType][i][2] += 2
                            ##self.ctrlCounts[seedType][i][3] += 1
                            if returnGenes:
                                rCtrlNames[seedType][i][1].extend([x for x in ctrlNames[i] if x.find(genomeToken) != -1])
                                type2 = [tuple(sorted(x)) for x in xpermutations.xuniqueCombinations(ctrlNames[i],2) if ''.join(x).find(genomeToken) != -1]
                                type3 = tuple(sorted(ctrlNames[i]))
                                rCtrlNames[seedType][i][2].extend(type2)
                                rCtrlNames[seedType][i][3].append(type3)

                if returnGenes:
                    for i in range(1,4):
                        assert len(rGeneNames[seedType][i]) == len(set(rGeneNames[seedType][i])),\
                               "ERROR: rGeneNames[%s] in miRNA(%s) has redundancy." % (i, self.name)
        if returnGenes:
            # store and return rGeneNames
            self.matchEvents = rGeneNames
            self.ctrlEvents  = rCtrlNames
            return (rGeneNames,rCtrlNames)
        
    def calcCtrlMeanStDv(self):
        """
        Iterate through the self.ctrlCount data and calculate the mean and stdDev for each
        seedType permutation set.
        """
        # Make sure we have already run self.countHitsInOrthos()
        assert self.matchCounts and self.ctrlCounts, \
               'ERROR:  It looks like we have not called countHitsInOrthos().'
        
        # Initialize data structure
        for seedType in _seedModels:
            self.ctrlMeanStd[seedType]  = [[None,None],[None,None],[None,None],[None,None]] # SEE self.ctrlStats def in __init__ for format definition
        
            
        # Iterate through self.ctrlCounts
        for seedType in self.ctrlCounts.keys():
            # Use numpy to record mean and StDv for each countCategory in self.ctrlCounts[seedType]
            for i in range(len(self.ctrlMeanStd[seedType])):
                meanStDv, mean = mathDefs.stdDv([x[i] for x in self.ctrlCounts[seedType]],kind='mean',df=1)
                self.ctrlMeanStd[seedType][i][0] = mean
                self.ctrlMeanStd[seedType][i][1] = meanStDv
                
    def calcCtrlMedianStDv(self):
        """
        Iterate through the self.ctrlCount data and calculate the MEDIAN and stdDev from the
        MEDIAN for each seedType permutation set.
        """
        # Make sure we have already run self.countHitsInOrthos()
        assert self.matchCounts and self.ctrlCounts, \
               'ERROR:  It looks like we have not called countHitsInOrthos().'
        
        # Initialize data structure
        for seedType in _seedModels:
            self.ctrlMedianStd[seedType]  = [[None,None],[None,None],[None,None],[None,None]] # SEE self.ctrlStats def in __init__ for format definition
        
        # Iterate through self.ctrlCounts
        for seedType in self.ctrlCounts:
            # Calc median and  medianStDv for each countCategory in self.ctrlCounts[seedType]
            for i in range(len(self.ctrlMedianStd[seedType])):
                medStDv, median = mathDefs.stdDv([x[i] for x in self.ctrlCounts[seedType]],kind='median',df=1)
                self.ctrlMedianStd[seedType][i][0] = median
                self.ctrlMedianStd[seedType][i][1] = medStDv
                
    def calcCtrlMedianMAD(self):
        """
        Iterate through the self.ctrlCount data and calculate the MEDIAN and the 
        Median Absolute Deviations for each seedType permutation set.
        """ 
        # Make sure we have already run self.countHitsInOrthos()
        assert self.matchCounts and self.ctrlCounts, \
               'ERROR:  It looks like we have not called countHitsInOrthos().'
        
        # Initialize data structure
        for seedType in _seedModels:
            self.ctrlMedianMAD[seedType]  = [[None,None],[None,None],[None,None],[None,None]] # SEE self.ctrlStats def in __init__ for format definition
        
        # Iterate through self.ctrlCounts
        for seedType in self.ctrlCounts:
            # Calc median and MAD for each countCategory in self.ctrlCounts[seedType]
            for i in range(len(self.ctrlMedianMAD[seedType])):
                MAD,median = mathDefs.medianAbsDev([x[i] for x in self.ctrlCounts[seedType]])
                self.ctrlMedianMAD[seedType][i][0] = median
                self.ctrlMedianMAD[seedType][i][1] = MAD
                
    def scoreSeedType(self, seedType, howManyOrthos, metric='medianStDv'):
        """
        scoreSeedType(self, seedType, howManyOrthos, metric='medianStDv')
        Returns the respective type of z-score.
        """
        metrics = {'meanStDv':self.ctrlMeanStd,
                   'medianStDv':self.ctrlMedianStd,
                   'medianMAD':self.ctrlMedianMAD,}
        
        
        assert seedType in _seedModels,\
               '\n\nERROR: Valid values of seedType include %s.' % (_seedModels)
        assert 0<=howManyOrthos<=3, \
               "\n\nERROR: howManyOrthos must be >= 0 and <= 3."
        assert metric in metrics,\
               "\n\nERROR: Valid values of metric include %s." % (metrics.keys())
        
        if not metrics[metric]:
            metMeths = {'meanStDv':self.calcCtrlMeanStDv(),
                        'medianStDv':self.calcCtrlMedianStDv(),
                        'medianMAD':self.calcCtrlMedianMAD(),}
            metMeths[metric]
        ##assert metrics[metric], \
               ##"\n\nERROR: %s has not been calculated yet." % (metric)
        
        real       = self.matchCounts[seedType][howManyOrthos]
        ##pprint(metrics[metric][seedType])
        ctrlCenter = metrics[metric][seedType][howManyOrthos][0]
        ctrlSpread = metrics[metric][seedType][howManyOrthos][1]
        
        return (float(real)-ctrlCenter)/ctrlSpread
        
    def reportGeneTargets(self,zScore=3,FDR=0.15):
        """
        Use supplied parameters to filter which seedTypes' "found" genes
        are returned.  Returns dict in form:
        keys(seedType):values(None,[geneNamesInOne,zScore,FDR],[geneNamesInTwo,zScore,FDR],[geneNamesInThree,zScore,FDR]).
        
        WARNING: as of now, this method COMPLETLY reRuns miRNA.countHitsInOrthos() with 
        returnGenes=True.  Do not run this over and over right now.
        """
        
        # initialize return Dict
        rDict = {}
        for seedType in _seedModels:
            rDict[seedType] = [None,[],[],[]]
        
        # Get geneNames
        geneNames = self.countHitsInOrthos(returnGenes=True)
        
        # Populate rDict
        for seedType in _seedModels:
            for orthoType in range(1,4):
                seedZScore = self.scoreSeedType(seedType,orthoType)
                seedFDR    = self.ctrlMedianStd[seedType][orthoType][0]/float(self.matchCounts[seedType][orthoType])
                if (seedZScore >= zScore) and (seedFDR < FDR):
                    rDict[seedType][orthoType].extend([geneNames[seedType][orthoType],seedZScore,seedFDR])
                else:
                    rDict[seedType][orthoType] = None
                
        return rDict
    
    def reportGeneTargetsFdrMedMeth(self,stdvLimit=2,consFdrThresh=0.25,divide=0):
        """
        Use supplied parameters to filter which seedTypes' "found" genes
        are returned.  Returns dict in form:
        keys(seedType):values(None,[geneNamesInOne,medianFDR,consFDR],[geneNamesInTwo,medianFDR,consFDR],[geneNamesInThree,medianFDR,consFDR]).
        
        **Eventually, this will allow the separate return of the grouped
        culicinae and Ag2Culic results based on the [divide] option.
        
        """
        
        if divide:
            pass # for now
        
        else:
            # initialize return Dict
            rDict = {}
            for seedType in _seedModels:
                rDict[seedType] = [None,[],[],[]]
            
            # Get geneNames
            if (self.matchEvents and self.ctrlEvents):
                geneNames = (self.matchEvents,self.ctrlEvents) # <- returns (geneNames,ctrlNames)
            else:
                geneNames = self.countHitsInOrthos(returnGenes=True) # <- returns (geneNames,ctrlNames)
            
            # Populate rDict
            for seedType in _seedModels:
                for orthoType in range(1,4):
                    realMatches    =  len(filterToken(geneNames[0][seedType][orthoType],'AGAP'))
                    ctrlMatches    =  [len(filterToken(x[orthoType],'AGAP')) for x in geneNames[1][seedType]]
                    
                    ctlrFDRstdv,ctlrFDRmed     =  calcFDRStats(realMatches,ctrlMatches) # <-returns(stdv,median)
                    
                    if (ctlrFDRstdv == None or ctlrFDRmed == None):
                        consSeedFDR = None
                    else:
                        consSeedFDR =  ctlrFDRmed+(stdvLimit*ctlrFDRstdv)
                    
                    if consSeedFDR == None:
                        rDict[seedType][orthoType] = None
                    else:
                        if round(consSeedFDR,2) <= consFdrThresh:
                            rDict[seedType][orthoType].extend([geneNames[0][seedType][orthoType],ctlrFDRmed,consSeedFDR])
                        else:
                            rDict[seedType][orthoType] = None
                    
            return rDict
    
    def purgeData(self):
        """
        Deletes self.matchData and self.ctrlData to reduce the size each miRNA obj.
        ONLY RUN AFTER matchData and ctrlData are no longer needed.
        """
        del self.matchData
        del self.ctrlData
            

#-#  Global defs  #-#
def calcFDRStats(realCount,ctrlCounts):

        FDRs = []
        
        for i in range(len(ctrlCounts)):
            ctrl = ctrlCounts[i]
            real = realCount
            
            if real == 0:
                return (None,None)
            else:
                fdr = float(ctrl)/real
                
            if fdr > 1:
                FDRs.append(1.0)
            else:
                FDRs.append(fdr)
                
        return mathDefs.stdDv(FDRs,'median')
        

        
def loadSeqs(fastaPathList):
    """
    Takes list of paths.  Returns single dict full of seqs found in the files.
    Converts softMasking to hard.
    """
    rDict = {}
    
    for path in fastaPathList:
        rDict.update(Fasta.load(path))
    
    bioDefs.softMaskDict2HardMask(rDict)
    return rDict

def loadOrthos(orthoPath, seqDict):
    """
    Takes path to ortho relation file.
    Returns list lof lists: [[orthoNameA, orthoNameB, orthoNameC]].
    Filters out any orthoLists that contain a gene that is missing from seqDict.
    """
    
    seqDictKeys = set(seqDict.keys())
    rList = []
    orthoPath = open(orthoPath, 'rU')
    originalCount = 0
    while 1:
        line = orthoPath.readline().strip()
        if not line: break
        originalCount+=1
        try:
            line = eval(line)
            assert type(line) == type([]), 'loadOrthos ERROR: type(line) != type([])'
            # Filter
            passed = True
            for gene in line:
                if gene not in seqDictKeys:
                    passed = False
                    break
            if passed:
                rList.append(line)
        except:
            rList.append(line.split('\t'))
    
    print '%s of %s original orthoRelations were left after filtering.' % (len(rList), originalCount) 
    return rList

def loadMiRNAs(miRNA_Path):
    """
    Takes fasta file of mature miRNAs.
    Returns dict.
    """
    
    return Fasta.load(miRNA_Path)

def filterOrthoSeqs(seqDict, orthoRelations):
    """
    Returns a dict with only those seqs that are named in orthoRelations.
    """
    flatList = sum([l for l in orthoRelations],[])
    
    rDict = {}
    
    for gene in flatList:
        rDict[gene] = seqDict[gene]
        
    return rDict

def getTokenFromNames(geneNames):
    """
    Returns the list of genome tokens from a list of gene names:
    ['AAEL','AGAP'] from geneNames=['AAEL0000567', 'AGAP0008908'].
    """
    toks = set()
    
    tokExp = re.compile('^\D+')
    for name in geneNames:
        match = tokExp.match(name)
        toks.add(match.group())
        
    return sorted(list(toks))
    
def filterToken(twoDList,token):
    rList =[]
    for i in range(len(twoDList)):
        for j in range(len(twoDList[i])):
            if twoDList[i][j].startswith(token):
                rList.append(twoDList[i][j])
    return list(set(rList))    
    


    
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
2009-08-15 -- Wrote miRNA.tallyHits()
2009-08-20 -- Wrote miRNA.countHitsInOrthos()
2009-08-20 -- Tested miRNA.tallyHits(): seems to work as expected
2009-08-20 -- Tested miRNA.countHitsInOrthos(): seems to work as expected
2009-08-26 -- Finished miRNA.loadSeqs().
2009-08-26 -- Finished miRNA.loadOrthos().
2009-08-27 -- Changed miRNA.tallyHits() -> initialization of self.ctrlData[seedType] from [set()]*number.  This creates many refs to ONE set().  NOT good.  Did it the hard way now.
2009-08-12 -- Added miRNA.reportGeneTargets(self,orthoType,zScore=3,FDR=0.15)
2009-09-14 -- Changed miRNA.countHitsInOrthos(self,returnGenes=False) --> added returnGenes=False to allow gene reporting.
"""

notes = """2009-05-25 -- Recording hits to each type of seed:
\t hit == dict{key=geneName : value=setOfModelsFound[model1,model2]}?
2009-05-25 -- Should move learnGenomes to MAIN so that I only do it once"""