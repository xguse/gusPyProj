from fjoin import FJoin


#========================= User Defined Variables =========================



# fjoin Paths
 
# codingBoundsOUT
path2fjIn1    = '/Users/biggus/Documents/MBGB/Rotations/James/Data/Sequence/Culex/culexProcessed/culexCodingBoundsOut.pos.txt '

# exom/codingExon file total
path2fjIn2    = '/Users/biggus/Documents/MBGB/Rotations/James/Data/Sequence/Culex/Culex_Exon_Coding_Reordered.txt'

# outfile
path2outFile  = '/Users/biggus/Documents/MBGB/Rotations/James/Data/Sequence/Culex/culexProcessed/culexBounds.pos.fjoin.txt'

# other fjoin arguments: [will need to change 5 & 6 when doing ]
""" set _columns1_ as correct variable from below """
col1_ValsUserDefined = '2,4,5'

col1_Vals4posOri_UpStrmReg = '2,9,10'   
col1_Vals4negOri_UpStrmReg = '2,11,12'
col1_Vals4posOri_DnStrmReg = ''   
col1_Vals4negOri_DnStrmReg = ''

columns1 =  col1_Vals4posOri_UpStrmReg # (for input 1) Modify numbers to represent HUMAN column numbers for chromosome, start, and end, respectively
columns2 = '2,6,7'  # (for input 1) See above
#==========================================================================


# Set up FJoin arguments and call FJoin
fjoinArgs = ["-1", path2fjIn1, 
             "-2", path2fjIn2, 
             "-s", "both", 
             "-o", path2outFile,
             "--columns1", columns1,
             "--columns2", columns2]

FJoin(fjoinArgs).go()
