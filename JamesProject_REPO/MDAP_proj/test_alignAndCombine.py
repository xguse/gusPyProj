from MDAP_defs import alignAndCombineMotifs
from TAMO import MotifTools

Motif = MotifTools.Motif



m = MotifTools.load('/Users/biggus/Documents/James/Collaborations/Campbell/data/Results_HyperGeoScreen/masked/Results_gGEMS/CCupAt4Days.6-8mers.gGEMS.top6.motifs.stdThresh.tmo')
w = [5.8952,
     5.6523,
     5.0585,
     4.9788,
     4.9678,
     4.7688]

toTmo = []
toTmo.append(alignAndCombineMotifs([m[0],m[1]],[w[0],w[1]]))
toTmo.append(alignAndCombineMotifs([m[0],m[4]],[w[0],w[4]]))
toTmo.append(alignAndCombineMotifs([m[1],m[4]],[w[1],w[4]]))
toTmo.append(alignAndCombineMotifs([m[2],m[3]],[w[2],w[3]]))

for e in toTmo:
    print e.oneletter

None