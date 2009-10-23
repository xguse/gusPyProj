from string import maketrans

#--------- Script Specific Function Definitions ---------------------


#--------------------------------------------------------------------



#========================= User Defined Variables =========================
inFile = '/Users/biggus/Documents/James/Data/2KB/2kb_Sequence/2kb_Culex/2kb_culexProcessed/culex2KBdownStreamTSS.softMasked.fas'
outFile = '/Users/biggus/Documents/James/Data/2KB/2kb_Sequence/2kb_Culex/2kb_culexProcessed/culex2KBdownStreamTSS.masked.fas'

#==========================================================================


inFile = open(inFile, 'rU').readlines()

transTab = maketrans('atcgn','NNNNN')

for i in range(0, len(inFile)):
    if inFile[i][0] == '>':
        pass
        header = inFile[i].split(':')
        inFile[i] = '>%s' % header[5]
    else:
        inFile[i] = inFile[i].translate(transTab)
        
outFile = open(outFile, 'w')

outFile.writelines(inFile)
        