

#========================= User Defined Variables =========================
inFile  = '/Users/biggus/Desktop/2KBup_collapsed_upstream_ex-conserved_mosquito-motifs_nr.te.hgp.ngt0.pVal.culled.txt'
outFile = '/Users/biggus/Desktop/2KBup_collapsed_upstream_ex-conserved_mosquito-motifs_nr.te.hgp.ngt0.pVal.culled.sorted.txt'

#==========================================================================

#--------- Script Specific Function Definitions ---------------------


#--------------------------------------------------------------------


inFile = open(inFile, 'rU').readlines()


for i in range (0, len(inFile)):
    line = inFile[i].split('\t')
    motifPair = [line[0], line[1]]
    motifPair.sort()
    line[0],line[1] = motifPair[0],motifPair[1]
    inFile[i] = '\t'.join(line)
    
    
outFile = open(outFile, 'w')
outFile.writelines(inFile)
outFile.close()
    