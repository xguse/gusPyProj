from TAMO.seq import Fasta


def groupMiRsBym2m8(miRNAs):
    """
    miRNAs = dict(k='miRname', v='miRseq')
    Returns seedDict = dict(k='m2m8Seq', v=[miRnames])
    """
    seedDict = {}
    
    for m in miRNAs:
        m2m8 = miRNAs[m][1:8]
        if m2m8 in seedDict:
            seedDict[m2m8].append(m)
        else:
            seedDict[m2m8] = [m]
        
    for each in seedDict:
        print '%s' % (', '.join(seedDict[each]))
        
    return seedDict


if __name__ == '__main__':
    miRNAs = '/Users/biggus/Documents/James/Data/Tu_miRNA/miRNAs/miRBase/mature.aga.fa'
    miRNAs = Fasta.load(miRNAs)
    sD = groupMiRsBym2m8(miRNAs)
    
