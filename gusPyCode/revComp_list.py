from gusPyCode.defs.JamesDefs import revComp


#========================= User Defined Variables =========================
inFile  = '/Users/biggus/Documents/James/Data/2KB/2kb_Sequence/2kb_Combo/2Kb_AllMosquitoes/MosqMotifs/MotifGroupPWMs/AllGroups/AllGroups_Aln.fwd.txt'
outFile = '/Users/biggus/Documents/James/Data/2KB/2kb_Sequence/2kb_Combo/2Kb_AllMosquitoes/MosqMotifs/MotifGroupPWMs/AllGroups/AllGroups_Aln.rvCmp.txt'

#==========================================================================

#--------- Script Specific Function Definitions ---------------------


#--------------------------------------------------------------------


inFile = map(lambda line : line.strip(), open(inFile, 'rU').readlines())


outData = []

for each in inFile:
    outData.append(revComp(each)+"\n")
    

outFile = open(outFile, 'w')
outFile.writelines(outData)

print "Done."