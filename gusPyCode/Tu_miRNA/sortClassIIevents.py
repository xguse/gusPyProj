from xpermutations import xuniqueCombinations


inFile  = '/Users/biggus/Documents/James/Data/Tu_miRNA/ResultsMiRNA_Targeting/2009_10_09run/2009_10_09runs.100bothCtrls.stats_1stdv.events.medFDRmeth.txt'
outFile = '/Users/biggus/Documents/James/Data/Tu_miRNA/ResultsMiRNA_Targeting/2009_10_09run/2009_10_09runs.stats_1stdv.medFDRmeth.sortedBy2s.new.txt'
genomeTokens = ['AAEL', 'AGAP', 'CPIJ'] 

inFile  = open(inFile, 'rU')
outFile = open(outFile, 'w')

def sortClass2sBySpcsPair(line,genomeTokens):
    """
    genomeTokens -> ['SpcA','SpcB','SpcC']
    
    writes results to outFile.
    """
    
    tokPairs = [tuple(sorted(x)) for x in xuniqueCombinations(genomeTokens,2)]    
    tokPairs.sort()
    
    data = {tokPairs[0]:[],
            tokPairs[1]:[],
            tokPairs[2]:[]}
    
    genePairs = eval(line[-1])
    
    for pair in genePairs:
        if   (tokPairs[0][0] in ''.join(pair)) and (tokPairs[0][1] in ''.join(pair)):
            data[tokPairs[0]].append(pair)
        elif (tokPairs[1][0] in ''.join(pair)) and (tokPairs[1][1] in ''.join(pair)):
            data[tokPairs[1]].append(pair)
        elif (tokPairs[2][0] in ''.join(pair)) and (tokPairs[2][1] in ''.join(pair)):
            data[tokPairs[2]].append(pair)
    
    line[-1:-1] = [str(len(data[tokPairs[0]])),str(len(data[tokPairs[1]])),str(len(data[tokPairs[2]]))]
    line[-1] = data[tokPairs[0]]+data[tokPairs[1]]+data[tokPairs[2]]
    # Write the counts to each pair       
    outFile.write('%s\t%s\n' % ('\t'.join(line[:-1]),line[-1]))
    

            






lines = []


# Organize by total or not
for line in inFile:
    if 'allPassedSeedsFor_2' in line:
        lines.append(line.strip('\n').replace("Seqs=", " : ").split(' : '))
    elif 'orthoType_2' in line:
        lines.append(line.strip('\n').replace("Seqs=", " : ").split(' : '))
    elif line.startswith('--'):
        pass
    
        
print 'Writing...'
for line in lines:
        sortClass2sBySpcsPair(line,genomeTokens)
    

    
print 'outFile: %s' % (outFile.name)
