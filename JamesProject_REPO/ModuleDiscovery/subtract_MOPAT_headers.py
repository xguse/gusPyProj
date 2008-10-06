from sets import Set



#========================= User Defined Variables =========================
inFile1 = '/Users/biggus/Documents/James/Data/2KB/2kb_Sequence/2kb_Combo/2Kb_AllMosquitoes/MosqMotifs/MotifGroupPWMs/AllGroups/ModuleData/AllGroups_MOPAT_int_VS_3hrDwnAll.headers.txt'
inFile2 = '/Users/biggus/Documents/James/Data/2KB/2kb_Sequence/2kb_Combo/2Kb_AllMosquitoes/MosqMotifs/MotifGroupPWMs/AllGroups/ModuleData/AllGroups_MOPAT_int_VS_3hrUpAll.headers.txt'

outFile = '/Users/biggus/Documents/James/Data/2KB/2kb_Sequence/2kb_Combo/2Kb_AllMosquitoes/MosqMotifs/MotifGroupPWMs/AllGroups/ModuleData/3hrDwnAll-3hrUpAll.txt'

#==========================================================================

#--------- Script Specific Function Definitions ---------------------


#--------------------------------------------------------------------

# Open both infiles and parse each by comma into lists at first
inFile1 = map(lambda line : line.strip().split(','), open(inFile1, 'rU').readlines())
inFile2 = map(lambda line : line.strip().split(','), open(inFile2, 'rU').readlines())

# convert each list's first index into a set after spliting on ';' for BOTH infiles
for i in range(0,len(inFile1)):
    index0_List = inFile1[i][0].split(';')
    inFile1[i][0] = Set(index0_List)
del i,index0_List
    
for i in range(0,len(inFile2)):
    index0_List = inFile2[i][0].split(';')
    inFile2[i][0] = Set(index0_List)

# Create lists to hold file2's sets
listOfSets2 = []
for each in inFile2:
    listOfSets2.append(each[0])
    
# var to hold subtracted sets
subtractedSets =[]

for each in inFile1:
    if each[0] not in listOfSets2:
        subtractedSets.append(str(each)+'\n')

print 'inFile1 had %i modules\ninfile2 had %i modules\noutFile will have %i modules' % (len(inFile1),len(inFile2),len(subtractedSets))        

outFile = open(outFile, 'w')
outFile.writelines(subtractedSets)

x=1