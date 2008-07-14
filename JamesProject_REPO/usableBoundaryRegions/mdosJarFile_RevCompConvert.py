from JamesDefs import revComp
import re

#--------- Script Specific Function Definitions ---------------------


#--------------------------------------------------------------------



#========================= User Defined Variables =========================
inFile =  '/Users/biggus/Documents/James/Data/2KB/2kb_Sequence/2kb_Combo/2kb_DmelMosquitoes/2KBup_CombinedDrosophilaAndCulexOrthologs.masked.fas'
outFile = '/Users/biggus/Documents/James/Data/2KB/2kb_Sequence/2kb_Combo/2kb_DmelMosquitoes/2KBup_CombinedDrosophilaAndCulexOrthologs.masked.rvCmp.fas'

#==========================================================================


inFile = map(lambda line : line.strip(), open(inFile, 'rU').readlines())

bp_10 = re.compile('^[ATGCN]{10,}')

for i in range(0,len(inFile)):
    if bp_10.search(inFile[i]):
        inFile[i] = "%snnn%s\n" % (inFile[i], revComp(inFile[i]))
    else:
        inFile[i] = inFile[i]+'\n'
        

        
    
outFile = open(outFile, 'w')

outFile.writelines(inFile)

print 'Done!'
        