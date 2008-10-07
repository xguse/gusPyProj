import JamesDefs
from defs_moduleByTableLookUp import makeFwdAndRevCompRegExObj_IUPAC
import re

#--------- Script Specific Function Definitions ---------------------


#--------------------------------------------------------------------



#========================= User Defined Variables =========================
inFile =  '/Users/biggus/Documents/James/Data/2KB/2kb_Sequence/2kb_Combo/2kb_DmelMosquitoes/Dmel-MosqMotifs/JAR_2KBupDmelAndMosqCombined_8mer.rvCmp.Reduced_gte3.motifs'
outFile = '/Users/biggus/Documents/James/Data/2KB/2kb_Sequence/2kb_Combo/2kb_DmelMosquitoes/Dmel-MosqMotifs/JAR_2KBupDmelAndMosqCombined_8mer.rvCmp.Reduced_gte3.nr.motifs'

#==========================================================================


inFile = map(lambda line : line.strip(), open(inFile, 'rU').readlines())


sortedMotifs = []

inFile_Len = len(inFile)
print "Starting number: %s" % inFile_Len


while inFile:
    
    # always inFile[0] bc i pop the tops later on so next motif always becomes inFile[0]
    motif = inFile[0].split('\t')
    motif = makeFwdAndRevCompRegExObj_IUPAC(motif[0])
    #print len(sortedMotifs)
    i = 0
    motifList = []
    
    while i < inFile_Len:
             
        if motif.search(inFile[i]):
            motifList.append(inFile.pop(i)+'\n')
            inFile_Len = len(inFile)
            i-=1
        
        i+=1
    sortedMotifs.append(motifList)

sortedMotifs_Len = len(sortedMotifs)

collapsedMotifList = []

for each in sortedMotifs:
    collapsedMotifList.append(each[0])

collapsedMotifList_Len = len(collapsedMotifList)
print "Final number: %s" % collapsedMotifList_Len

x=1

outFile = open(outFile, 'w')
outFile.writelines(collapsedMotifList)

X=1
