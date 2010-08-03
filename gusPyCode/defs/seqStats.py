import sys
from TAMO.seq import Fasta





def calcStats(fastaPath):
    seqFile = Fasta.load(fastaPath)
    combinedSeq = ''
    
    for each in seqFile:
        combinedSeq += seqFile[each]
    
    combinedSeq= combinedSeq.upper()
    
    seqs       = len(seqFile)
    totNucs    = len(combinedSeq)
    aCnt       = combinedSeq.count('A')
    cCnt       = combinedSeq.count('C')
    gCnt       = combinedSeq.count('G')
    tCnt       = combinedSeq.count('T')
    nCnt       = combinedSeq.count('N')
    nonNs      = aCnt+cCnt+gCnt+tCnt
    n2tot      = float(nCnt)/len(combinedSeq)
    n2nonN     = float(nCnt)/nonNs
    percentGC  = (float(gCnt)+cCnt)/nonNs
    
    
    
    return {'seqLen':seqs,
            'totNucs':totNucs,
            'aCnt':aCnt,
            'cCnt':cCnt,
            'gCnt':gCnt,
            'tCnt':tCnt,
            'nCnt':nCnt,
            'nonNs':nonNs,
            'n2tot':n2tot,
            'n2nonN':n2nonN,
            'percentGC':percentGC}

    
    
if __name__ == '__main__':
    assert len(sys.argv[1:]) == 1 , "usage: %s [fastaFile]" % (sys.argv[0].split('/')[-1])


    
    print '\nfileName supplied is %s' % (sys.argv[1])
    
    stats = calcStats(sys.argv[1])
    
    print '''Total Sequences:\t\t%s
Total Nucleotides:\t\t%s
Total Non-N Nucleotides:\t%s
Total N Nucleotieds:\t%s
Ns/tot:\t\t\t%s
Ns/non-Ns:\t\t\t%s
A count:\t\t\t%s
C count:\t\t\t%s
G count:\t\t\t%s
T count:\t\t\t%s
Percent GC:\t\t\t%s''' % (stats['seqLen'],stats['totNucs'],stats['nonNs'],stats['nCnt'],stats['n2tot'],stats['n2nonN'],stats['aCnt'],stats['cCnt'],stats['gCnt'],stats['tCnt'],stats['percentGC'])
    
    
    
    
