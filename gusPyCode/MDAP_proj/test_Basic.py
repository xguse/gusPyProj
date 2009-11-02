from gusPyCode.MDAP_proj.MD_wrappers import TamoWrap

# TamoWrap -> def __init__(self, optionsObj, posArgs):

optionsObj = {'kmerSize':7,'kmerRange':'6,10'} # need to change the use of this to be how mortals think of numbers
posArgs    = ['/Users/biggus/Documents/James/Data/2KB/2kb_Sequence/2kb_Anopheles/2KBupTSS_goodAffyAGAPsFastasOUT.masked.nr.seqsShuffled.fas','/Users/biggus/Documents/James/Data/ReClustering/kmedsPear33Clus50x_2/Clus2_247genes.genes.txt']
outFile    = '/Users/biggus/Documents/James/Data/ReClustering/kmedsPear33Clus50x_2/Clus2_247genes.6-8mers.txt'

tWrap = TamoWrap(optionsObj,posArgs)
tWrap.go()

#for l in tWrap.toFile:
    #print l.strip()

print 'num of keepers = %s' % (len(tWrap.output))

outFile = open(outFile, 'w')
outFile.writelines(tWrap.toFile)