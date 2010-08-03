from gusPyCode.defs import motifs

xmsPath = '/Users/biggus/Documents/James/Writings_Talks/Manuscripts/AaTxDelta_RNAseq/creDiscovery/NTO/LSvLB.nuls.Tx.2000.xms'
xmsMotifs = motifs.parseXMSfile(xmsPath)
motifs.writeXMSfile(xmsMotifs,'/Users/biggus/Documents/James/Writings_Talks/Manuscripts/AaTxDelta_RNAseq/creDiscovery/NTO/LSvLB.nuls.Tx.2000.2.xms')

None