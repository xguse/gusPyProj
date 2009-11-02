from time import time
from TAMO.MotifTools import Motif
from TAMO.seq import Fasta
from gusPyCode.defs.statsDefs import cumHypergeoP



class SeqMap:
    """
    Class to create/store hit maps of a dictionary of TAMO motif objects.
    Default score threshold is set to 0.5 to give the user flexibility downstream.
    Constructor: SeqMap(seqName, sequence, motifDict,thresh=0.5)
    """
    def __init__(self, seqName, sequence, motifDict,thresh=0.5,window=200):
        
        self.name       = seqName
        self.motifSets  = {}        # key<(winSize,threshUsed)>, Value<FrozenSet-of-FrozenSets>
        self.motifPairs = {}        # key<(winSize,threshUsed)>, Value<FrozenSet-of-FrozenSets>
        self.locMap     = []
        self.motifs     = {}
        
        for m in motifDict:
            assert motifDict[m].__class__ == Motif().__class__, \
                   'motifDict contains at least one non-TAMO motif (%s)' % (m)
            
            maxScore = motifDict[m].maxscore
            hits     = motifDict[m].scan(sequence,factor=thresh)
            
            self._updateMap(motifDict[m],m,hits)
            
            
        self.locMap.sort()
        #if   window: self.buildMotifPairs(window)
        ##elif window >  2: self._buildMotifSets(window)
        ##elif window <  2: print 'ALERT: SeqMap argument(window) set below 2. No combo sets produced.'
        
            
    def _updateMap(self, motif, motifID, motifScan):
        """
        Converts motifScan tuple into individual locations for each hit.
        Location format: [start, stop, motifID, str(motif), score, matrixMaxScore, hitSeq]
        """
        hitSet = []
        
        for i in range(len(motifScan[0])):
            
            hitSet.append([motifScan[1][i][0],  # start
                           motifScan[1][i][1],  # stop
                           motifID,             # name
                           str(motif),          # one letter representation of motif
                           motifScan[2][i],     # score hit got
                           motif.maxscore,      # max score posible from motif matrix
                           motifScan[0][i]])    # hitSeq
        
        self.locMap.extend(hitSet)
        
        # Add motif to motifDict if any hits occur.
        if hitSet:
            self.motifs[motifID]=motifScan
        #print 'Updated seqMap(%s) for %s' % (self.name,motifID)
    
    def buildMotifPairs(self, window=200,thresh=0.75):
        """
        Constructs a FrozenSet-of-FrozenSets representing all sets of non-redundant pairs of motifIDs found
        within the specified window for a SeqMap obj.
        
        Motif starting positions are not allowed to be less than 5 bp apart.
        """
        # //// Test Results \\\\
        # _x_ Only Seeds windows with locs of right thresh
        # _x_ Only populates windows with locs of right thresh
        # _x_ Correctly excludes overlapping motifPairs at right length
        # _x_ Double motifs are recorded correctly.
        
        mapPairSet = set()
        
        # Generate the Set of motifID pairs representing all pairs of non-redundant motifIDs found within the specified window for a SeqMap obj
        for i in range(len(self.locMap)):
            # Enforce thresh level
            threshTest1 = float(self.locMap[i][4])/self.locMap[i][4]
            if threshTest1 < thresh:
                continue
            winLocs = []
            end     = self.locMap[i][0] + window
            # Generate each window's List of motif locations
            for jLoc in self.locMap[i:len(self.locMap)]:
                threshTest2 = float(jLoc[4])/jLoc[5]
                if threshTest2 < thresh: continue  ## Only include locations that meet thresh level
                if jLoc[0] >    end: break     ## Stop populating window at end of window len
                winLocs.append(jLoc)
            if winLocs:
                # Extract pairs of motifIDs with starting pos at least 5 bp apart.
                winPairs = []
                primeMotif = winLocs[0]
                for i in range(1,len(winLocs)):
                    if winLocs[i][0] - primeMotif[0] >= 5:   # This is bc I am INCLUDING the 1st motif's 1st pos.
                        if winLocs[i][2] == primeMotif[2]:
                            winPairs.append(frozenset([primeMotif[2],'double']))  # insures that double motifs will be distinuishable from possible erroneously included single motifs
                            continue
                        winPairs.append(frozenset([primeMotif[2],winLocs[i][2]]))
                
                if winPairs: mapPairSet.update(winPairs)
        self.motifPairs[(window,thresh)] = frozenset(mapPairSet)

    
    def _buildMotifSets(self, windowSize):
        """
        Constructs a FrozenSet-of-FrozenSets representing all sets of non-redundant motifIDs found
        within the specified window for a SeqMap obj.
        
        THIS NOT IMPLEMENTED YET -> Motif starting positions are not allowed to be less than 4 bp apart.
        """

        '''!!  May want to only keep best scoring motif of those with same oneletter string  !!'''
        
        mapMotifSet = set()
        
        # Generate the Set of motifID Sets representing all sets of non-redundant motifIDs found within the specified window for a SeqMap obj
        for i in range(len(self.locMap)):
            winLocs = []
            end     = self.locMap[i][0] + windowSize
            # Generate each window's List of motif locations
            for jLoc in self.locMap[i:len(self.locMap)]:
                if jLoc[0] > end:
                    break
               # if 
                winLocs.append(jLoc)
            # Convert window location list to window motifID Set.
    ## !! _May_ want to only keep best scoring motif of those with same oneletter string  !!
            winMotifSet = set()
            for each in winLocs:
                winMotifSet.add(each[2])
            # Add this window set to the mapSet
            mapMotifSet.add(frozenset(winMotifSet))
            
        # Add this mapSet to the dict of mapSets with it's winSize as key
        self.motifSets[windowSize] = frozenset(mapMotifSet)
            
            
class MapLib:
    """
    Class to batch-instantiate, store, and work with multiple SeqMap objects.
    Constructor: MapLib(fastaSeqs,motifDict,thresh=0.5)
    """
    def __init__(self,fastaSeqs, motifDict, thresh=0.5,window=200):
        
        self.seqMaps = {}
        
        # Get seqs from fasta
        assert type(fastaSeqs) == type('string') or type(fastaSeqs) == type({}),\
               'MapLib arg(fastaSeqs) must be string pointing to file or a seqDict.'
        if type(fastaSeqs) == type('string'):
            seqs = Fasta.load(fastaSeqs)
        elif type(fastaSeqs) == type({}):
            seqs = fastaSeqs
        
        # Instantiate a SeqMap obj for each seq in seqs
        c = 0
        for k in seqs:
            c += 1
            assert c <= 250
            realT1 = time()
            self.seqMaps[k] = SeqMap(k, seqs[k], motifDict, thresh=thresh, window=window)
            realT2 = time()
            print '%.4f\t%s' % (realT2-realT1,c)
    
    def getMotifCombos4GeneList(self, geneList, comboSize=2):
        pass
    
    def getMotifPairs4GeneList(self, geneList, window=200, thresh=0.75):
        """
        Returns non-redundant set of motif pairs found in the supplied list of geneNames.
        """
        # //// To Address \\\\
        # ___ Warn if geneName in geneList but NOT in mapLib
        
        
        motifSets = set()
        
        # Check to see if the data requested exists in each SeqMap obj; if not, produce it. 
        for sMap in self.seqMaps:
            keys = self.seqMaps[sMap].motifPairs.keys()
            if (window,thresh) in keys:
                continue
            else:
                self.seqMaps[sMap].buildMotifPairs(window,thresh)
        
        for sMap in self.seqMaps:
            if self.seqMaps[sMap].name in geneList:
                motifSets.update(self.seqMaps[sMap].motifPairs[(window,thresh)])
        
        return motifSets
    
    def calcPairsCumHG_pVal(self, geneList, motifPair, window=200, thresh=0.75):
        """
        Returns the cumulative hypergeometric p-value for a motifPair in the geneList compared
        to the background genes.
        """
        
        n = 0  # numb of positives in population
        i = 0  # numb of positives in sample
        m = 0  # numb of negatives in population
        N = 0  # sample size
        
        
        for sMap in self.seqMaps:
            if motifPair in self.seqMaps[sMap].motifPairs[(window,thresh)]:
                n += 1
                if self.seqMaps[sMap].name in geneList:
                    i += 1
            else:
                m += 1
        N = len(geneList)
        
        return {'pair':motifPair,'n':n,'i':i,'m':m,'N':N,'p':cumHypergeoP(n,i,m,N)}
    
