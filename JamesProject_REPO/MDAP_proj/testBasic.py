from MD_wrappers import TamoWrap

# TamoWrap -> def __init__(self, optionsObj, posArgs):

optionsObj = {'kmerSize':9,'kmerRange':None}
posArgs    = ['/Users/biggus/Documents/James/Data/2KB/2kb_Sequence/2kb_Anopheles/2KBupTSS_goodAffyAGAPsFastasOUT.masked.nr.fas','/Users/biggus/Documents/James/Data/ReClustering/kMedsPear33Clus50x/kMedsPear33Clus50x-1.genes.txt']


tWrap = TamoWrap(optionsObj,posArgs)
tWrap.go()

for l in tWrap.toFile:
    print l.strip()

print 'num of keepers = %s' % (len(tWrap.output))
