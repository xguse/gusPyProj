
"""
Takes file in form:
AAEL000001-RA:	AAEL000001-RA	AAEL007286-RA
AAEL000002-RA:	AAEL000002-RA	AAEL004300-RA

Returns file filtered for those that have at least [atLeastNum] Tx matches excluding self.
[onlyGenes=True] option in loadLines(filePath,onlyGenes=None) treats all transcripts as a
single entity ~Gene.

"""
# <userSettings>
inFile     = ''
outfile    = ''
atLeastNum = 10
onlyGenes  = None
# </userSettings>


# <defs>
def loadLines(filePath,onlyGenes=None):
    """
    Takes filePath. Returns Dict with probeName as Key and list of Matches as Value.
    """
    
    matchSets = {}
    for line in open(filePath, 'rU'):
        line    = line.strip().split()
        probe   = line[0].replace(':','')
        matches = line[1:]
        if onlyGenes:
            # Remove '-R*' from the Tx Name (AGAP000000-RA) and reduce to nr
            # list of geneNames rather than TxNames.
            for i in range(len(matches)):
                matches[i] = matches[i][:-3]
            matches = list(set(matches))
            matches.sort()
                    
        matchSets[probe] = matches
    return matchSets
#</defs>

# # # # # # # # --Main-- # # # # # # # #

matchSets = loadLines(inFile,onlyGenes=onlyGenes)



