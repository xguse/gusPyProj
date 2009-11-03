from gusPyCode.MDAP_proj.MDAP_defs import findBestPairAlignments
from TAMO.MotifTools import Motif,load


iFiles = ['/Users/biggus/Documents/James/Collaborations/Campbell/data/Results_HyperGeoScreen/masked/Results_gGEMS/CCupAt4Days.6-8mers.gGEMS.top6.motifs.stdThresh.tmo']
motifs = []

for i in range(len(iFiles)):
    motifs.extend(load(iFiles[i]))


    
mat = findBestPairAlignments(motifs, 1)

None
