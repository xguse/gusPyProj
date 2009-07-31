
histData = '/Users/biggus/Documents/James/Data/Tu_miRNA/MDOSoutPut/MDOS_kmer_IDs/AaAgAndAgCq_7mers_vs_AgMiRNAs.allMotifsNoAmbig.histData.txt'
histData = map(lambda l: l.strip().split('\t'), open(histData,'rU').readlines())

oFile = '/Users/biggus/Documents/James/Data/Tu_miRNA/MDOSoutPut/MDOS_kmer_IDs/AaAgAndAgCq_7mers_vs_AgMiRNAs.allMotifsNoAmbig.histData.trpsd.txt'
oFile = open(oFile, 'w')

for i in range(len(histData[0])):
    line = ''
    try:
        if histData[0][i]:
            line+='%s\t' % (histData[0][i])
    except:
        line+='\t'
    try:
        if histData[1][i]:
            line+='%s\n' % (histData[1][i])
    except:
        line+='\n'
    
    print line[:-1]
    oFile.write(line)