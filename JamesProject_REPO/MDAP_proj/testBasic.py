from MD_wrappers import TamoWrap

# TamoWrap -> def __init__(self, optionsObj, posArgs):

optionsObj = {'kmerSize':7,'kmerRange':None}
posArgs    = ['/Users/biggus/Documents/James/Data/2KB/2kb_Sequence/2kb_Anopheles/2KBupTSS_goodAffyAGAPsFastasOUT.masked.nr.fas','/Users/biggus/Documents/James/Data/ReClustering/kmedsPear33Clus50x_2/kmedsPear33Clus50x_2-2.genes.txt']
outFile    = '/Users/biggus/Documents/James/Writings_Talks/Grants/09_Feb/PrelimData_Grant_Feb09/Clus2_7merSearch-0.01.txt'

tWrap = TamoWrap(optionsObj,posArgs)
tWrap.go()

#for l in tWrap.toFile:
    #print l.strip()

print 'num of keepers = %s' % (len(tWrap.output))

outFile = open(outFile, 'w')
outFile.writelines(tWrap.toFile)