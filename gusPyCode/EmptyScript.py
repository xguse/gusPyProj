from doug_hypergeometric import hyperGeoPvalue
from gusPyCode.defs import JamesDefs




#--------- Script Specific Function Definitions ---------------------


#--------------------------------------------------------------------



#========================= User Defined Variables =========================
inFile = '/Users/biggus/Documents/MBGB/Rotations/James/Data/mdosJAR_testing/JAR_2KBupAedesAnopheles_7mer.rvCmp.smt2.sortedMotifsOnly.motifs'
#outFile = '/Users/biggus/Documents/MBGB/Rotations/James/Data/mdosJAR_testing/JAR_2KBupAedesAnopheles_7mer.rvCmp.smt2.sorted.motifs'

#==========================================================================


inFile = map(lambda line : line.strip(), open(inFile, 'rU').readlines())

nrList = JamesDefs.nrListBySets(inFile)




        
        
x=1
        
        