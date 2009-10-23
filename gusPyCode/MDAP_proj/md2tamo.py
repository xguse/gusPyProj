from gusPyCode.MDAP_proj.MDAP_defs import loadMotifsFromOutFile
from TAMO import MotifTools


mdOutFiles = ['/Users/biggus/Documents/James/Data/ReClustering/PrelimData_Grant_Feb09/RandSplitFastas/AceResults/Clus2_247gene_0.8_Apr16_14-46-33.ace.1.txt',
             '/Users/biggus/Documents/James/Data/ReClustering/PrelimData_Grant_Feb09/RandSplitFastas/AceResults/Clus2_247gene_0.8_Apr16_14-46-33.ace.2.txt',
             '/Users/biggus/Documents/James/Data/ReClustering/PrelimData_Grant_Feb09/RandSplitFastas/AceResults/Clus2_247gene_0.8_Apr16_14-46-33.ace.3.txt']


for mdFile in mdOutFiles:
    motifs = loadMotifsFromOutFile(mdFile,'list') # ['Meme', 'AlignAce', 'MDscan', 'Weeder','list']
    MotifTools.save_motifs(motifs,mdFile+'.tmo')
    
print 'Done.'

