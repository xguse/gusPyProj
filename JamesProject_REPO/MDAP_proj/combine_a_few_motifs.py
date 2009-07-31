from TAMO.MotifTools import load,sum
import MDAP_defs as md


motifs = load('/Users/biggus/Documents/James/Collaborations/Campbell/data/Results_HyperGeoScreen/masked/Results_gGEMS/CCupAt4Days.6-8mers.gGEMS.top6.motifs.stdThresh.tmo')

# Weights __MUST__ be in same order as motifs are loaded.
weights = [1,
           2,
           3,
           4,
           -1,
           5]

assert len(motifs) == len(weights), \
       'Error: len(motifs) MUST equal len(weights).'

motifs = zip(*[motifs, weights])
motifs.sort(key=lambda x: x[1])
motifs = zip(*motifs)

aligned = md.alignSimilarMotifs(motifs[0])
combined = sum(aligned,zip(*motifs[1]))

None

