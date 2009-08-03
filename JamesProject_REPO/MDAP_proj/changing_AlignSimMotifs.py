from MDAP_defs import findBestPairAlignments
from MDAP_defs import alignAndCombineMotifs
from TAMO import MotifTools
Motif = MotifTools.Motif

def alignSimilarMotifs(similarMotifs, minoverlap=6):
    """
    Takes list of similar motifs. List should be sorted by weight or pvalue since this def
    basically aligns each motif to the top motif.  Its kind of a progressive pairwise alignment.
    """
    alignedMotifs = []
    matrix = findBestPairAlignments(similarMotifs, minoverlap=minoverlap, verbose=None)
    
    # Find longest neg-offset to motif0
    lNegOff = 0
    for i in range(1,len(similarMotifs)):
        if matrix[0][i][3] < lNegOff:
            lNegOff = matrix[0][i][3]
    # Adjust left padding of motif0 for longest negOffset
    alignedMotifs.append(similarMotifs[0][lNegOff,similarMotifs[0].width])
    None
    
    # Adjust orientation and left padding for each motif
    for i in range(1,len(similarMotifs)):
        if matrix[0][i][3] < 0:                                   # If neg offset
            rcMotif = similarMotifs[i].revcomp()                  #   revComp motif_i 
            lPad = lNegOff-matrix[0][i][3]                        #   remember -> neg - neg = closer to 0
            alignedMotifs.append(rcMotif[lPad,rcMotif.width])
        elif matrix[0][i][3] > 0:                                 # If pos offset
            lPad = lNegOff-matrix[0][i][3]                        #   remember -> neg - pos = farther from 0
            alignedMotifs.append(similarMotifs[i][lPad,similarMotifs[i].width])
        elif matrix[0][i][3] == 0:                                # If no offset
            alignedMotifs.append(similarMotifs[i][lNegOff,similarMotifs[i].width])
    
    # Add right padding to match longest motif from above
    # Find Longest motif
    lMotifLen = 0
    for mtf in alignedMotifs:
        if mtf.width > lMotifLen:
            lMotifLen = mtf.width
    for i in range(len(alignedMotifs)):
        if alignedMotifs[i].width < lMotifLen:
            alignedMotifs[i] = alignedMotifs[i][0,lMotifLen]
        
        
    
    return alignedMotifs


m = MotifTools.load('/Users/biggus/Documents/James/Collaborations/Campbell/data/Results_HyperGeoScreen/masked/Results_gGEMS/CCupAt4Days.6-8mers.gGEMS.top6.motifs.stdThresh.tmo')
w = [5.8952,
     5.6523,
     5.0585,
     4.9788,
     4.9678,
     4.7688]

twoFive = [[m[0],m[1],m[4]],[w[0],w[1],w[4]]]

alndMotifs = alignSimilarMotifs(twoFive[0], minoverlap=4)
for m in alndMotifs:
    print m.oneletter
    
sumdMotif = MotifTools.sum(alndMotifs)

#bKmers = sumdMotif.bogus_kmers()
#for k in bKmers:
    #print k
None