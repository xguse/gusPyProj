
"""
Takes file in form:
AAEL000001-RA:	AAEL000001-RA	AAEL007286-RA
AAEL000002-RA:	AAEL000002-RA	AAEL004300-RA

Returns file filtered for those that have at least [atLeastNum] Tx matches excluding self.
[onlyGenes=True] option in loadLines(filePath,onlyGenes=None) treats all transcripts as a
single entity ~Gene.

"""
# <userSettings>
inFile     = '/Users/biggus/Documents/James/Data/MicroArrayPrep/Aedes/Probes/ProbeBlastResults/Aedes.agilent.probes.nr.blastn_low.AedesTxs.30sets.txt'
outfile    = '/Users/biggus/Documents/James/Data/MicroArrayPrep/Aedes/Probes/ProbeBlastResults/Aedes.agilent.probes.nr.blastn_low.AedesTxs.30sets.atLeast10.toTheGene.txt'
atLeastNum = 10
onlyGenes  = True
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

rtnProbes = []

for probe in matchSets:
    if onlyGenes:
        assert probe[:-3] in matchSets[probe], \
               "ERROR: Probe %s does not seem to have matched itself." % (probe)
    else:
        assert probe in matchSets[probe], \
               "ERROR: Probe %s does not seem to have matched itself." % (probe)
    
    if len(matchSets[probe]) >= atLeastNum+1: # we will not count its match to itself
        rtnProbes.append(probe)

rtnProbes.sort()


# write out filtered List:
outfile = open(outfile, 'w')
for probe in rtnProbes:
    outfile.write('%s:\t%s\n' % (probe,'\t'.join(matchSets[probe])))
    
print 'Written %s sets with at least %s matches to the probe.' % (len(rtnProbes), atLeastNum)

