

#--------- Script Specific Function Definitions ---------------------


#--------------------------------------------------------------------



#========================= User Defined Variables =========================
convertionKey = '/Users/biggus/Documents/James/Data/2KB/2kb_Sequence/2kb_Drosophila/FBgn_2_CG.nr.txt'
fastaFile     = '/Users/biggus/Documents/James/Data/2KB/2kb_Sequence/2kb_Drosophila/drosophila2KBupStreamTSS.masked.fas'
outFile       = '/Users/biggus/Documents/James/Data/2KB/2kb_Sequence/2kb_Drosophila/drosophila2KBupStreamTSS.nameConvertion.masked.fas'

#==========================================================================


convertionKey = map(lambda line: line.strip(), open(convertionKey, 'rU').readlines())
fastaFile     = open(fastaFile, 'rU').readlines()



for i in range(0, len(convertionKey)):
    convertionKey[i] = convertionKey[i].split('\t')
    
cvrsnDict = {}

for each in convertionKey:
    cvrsnDict[each[0]] = each[1]
    
#print convertionKey[0:3]
#print fastaFile[0:6]    
x=0


for i in range(0, len(fastaFile)):
    if fastaFile[i][0] == '>':
        old = fastaFile[i][1:12]
        fastaFile[i] = fastaFile[i].replace(old, "%s\t%s" % (cvrsnDict[old],old))
        print 'old:%s new:%s' % (old, fastaFile[i][1:12])
        

x = 1            

        
outFile = open(outFile, 'w')

outFile.writelines(fastaFile)
        