from gusPyCode.defs.doug_hypergeometric import hyperGeoPvalue
from gusPyCode.defs import JamesDefs




#--------- Script Specific Function Definitions ---------------------


#--------------------------------------------------------------------



#========================= User Defined Variables =========================
inFile = '/Users/biggus/Documents/James/Data/OrthologDefs/nrOrthos/CqMosqOrthos_count3s.txt'
#outFile = '/Users/biggus/Documents/MBGB/Rotations/James/Data/mdosJAR_testing/JAR_2KBupAedesAnopheles_7mer.rvCmp.smt2.sorted.motifs'

#==========================================================================


inFile = map(lambda line : line.strip(), open(inFile, 'rU').readlines())

c = 0
w = 0
for i in range(0,len(inFile)):
    print 'i = %i' % (i)
    testList = inFile[i:i+3]
    if len(testList) == 3: 
        if testList[0] == testList[1] == testList[2]:
            c+=1
    else:
        w+=1
print "-----\nc = %i\nw = %i" % (c,w)
    


        
        
x=1
        
        