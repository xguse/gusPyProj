from gusPyCode.defs import JamesDefs
import re

#--------- Script Specific Function Definitions ---------------------
def makeFwdAndRevCompRegExObj_IUPAC(motif, equals=0):
    import re
    
    motifPair = [motif, JamesDefs.revComp(motif)]
        
    targetContainsMotif  = '(%s|%s)' % (motifPair[0], motifPair[1])
    targetISMotif        = '^(%s|%s)&' % (motifPair[0], motifPair[1])
    
    if equals == 1:
        motif = targetISMotif
    elif equals == 0:
        motif = targetContainsMotif
    
    fwdRevComp_regExObj = re.compile(motif, re.IGNORECASE)
    return fwdRevComp_regExObj


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

#x=1

outFile = open(outFile, 'w')
outFile.writelines(collapsedMotifList)

#X=1


def reduceRevCompMotifs(listOfMotifs_IUPAC):
    reducedMotifs = []
    
    pass