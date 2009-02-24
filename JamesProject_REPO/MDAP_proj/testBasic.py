from MD_wrappers import TamoWrap

# TamoWrap -> def __init__(self, optionsObj, posArgs):

optionsObj = {'kmerSize':7,'kmerRange':'6,9'} # need to change the use of this to be how mortals think of numbers
posArgs    = ['/Users/biggus/Documents/James/Data/2KB/2kb_Sequence/2kb_Anopheles/2KBupTSS_goodAffyAGAPsFastasOUT.masked.nr.fas','/Users/biggus/Documents/James/Data/DoubleBloodMeal/DBMclusterGenes.txt']
outFile    = '/Users/biggus/Documents/James/Data/DoubleBloodMeal/DBMgenesKmerAnalysis.6-8mers.txt'

tWrap = TamoWrap(optionsObj,posArgs)
tWrap.go()

#for l in tWrap.toFile:
    #print l.strip()

print 'num of keepers = %s' % (len(tWrap.output))

outFile = open(outFile, 'w')
outFile.writelines(tWrap.toFile)