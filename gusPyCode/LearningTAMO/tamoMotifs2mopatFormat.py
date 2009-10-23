import cPickle
from TAMO.MotifTools import load as loadTMOs
from gusPyCode.defs.bioDefs import convertTAMOs2MOPATs

outFile  = '/Users/biggus/Documents/James/Data/ReClustering/PrelimData_Grant_Feb09/Clus2_247genes.Ace_gGEMS_MEME.mopat.motifs.txt'
tPkls    = None #['',]
tmoFiles = ['/Users/biggus/Documents/James/Data/ReClustering/PrelimData_Grant_Feb09/RandSplitFastas/AceResults/Clus2_247gene_0.8_Apr16_14-46-36.ace.2.txt.tmo',
            '/Users/biggus/Documents/James/Data/ReClustering/PrelimData_Grant_Feb09/RandSplitFastas/MemeResults/Clus2_247gene_0.8_Apr16_14-46-36.meme.txt.tmo',
            '/Users/biggus/Documents/James/Data/ReClustering/PrelimData_Grant_Feb09/RandSplitFastas/MemeResults/Clus2_247gene_0.8_Apr16_14-46-33.meme.txt.tmo',
            '/Users/biggus/Documents/James/Data/ReClustering/PrelimData_Grant_Feb09/RandSplitFastas/AceResults/Clus2_247gene_0.8_Apr16_14-46-36.ace.3.txt.tmo',
            '/Users/biggus/Documents/James/Data/ReClustering/PrelimData_Grant_Feb09/RandSplitFastas/AceResults/Clus2_247gene_0.8_Apr16_14-46-36.ace.1.txt.tmo',
            '/Users/biggus/Documents/James/Data/ReClustering/PrelimData_Grant_Feb09/RandSplitFastas/AceResults/Clus2_247gene_0.8_Apr16_14-46-33.ace.3.txt.tmo',
            '/Users/biggus/Documents/James/Data/ReClustering/PrelimData_Grant_Feb09/RandSplitFastas/AceResults/Clus2_247gene_0.8_Apr16_14-46-33.ace.2.txt.tmo',
            '/Users/biggus/Documents/James/Data/ReClustering/PrelimData_Grant_Feb09/RandSplitFastas/AceResults/Clus2_247gene_0.8_Apr16_14-46-33.ace.1.txt.tmo',
            '/Users/biggus/Documents/James/Data/ReClustering/PrelimData_Grant_Feb09/Clus2_247genes.6-8mers.gGEMS.tmo']

motifs = []
if tPkls != None: 
    for p in tPkls:
        motifs.extend(cPickle.load(open(p, 'r')))

if tmoFiles != None:
    for t in tmoFiles:
        motifs.extend(loadTMOs(t))

motifs = convertTAMOs2MOPATs(motifs)
    
# print to file
keys = sorted(motifs.keys(),key=lambda x: int(x.split('_')[0]))

outFile = open(outFile, 'w')
for k in keys:
    toOut = '>%s\tBlank\n' % (k)
    for i in range(len(motifs[k])):
        toOut += '%s\n' % (' '.join(map(lambda x: str(x),motifs[k][i])))
    outFile.write(toOut)
    
outFile.close()

print 'Done.'    

    
    