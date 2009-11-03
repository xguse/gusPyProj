from gusPyCode.MDOSX_proj.MDOSX_defs.defs import combineOrthoTabs
from pprint import pprint
from time import time

#file1 = '/Users/biggus/Documents/James/Data/OrthologDefs/nrOrthos/!nrOrthoDefs/Aaeg_Agam_1-to-1.nr.txt'
#file2 = '/Users/biggus/Documents/James/Data/OrthologDefs/nrOrthos/!nrOrthoDefs/Aedes_Culex_1-to-1.nr.txt'
#file3 = '/Users/biggus/Documents/James/Data/OrthologDefs/nrOrthos/!nrOrthoDefs/Culex_Agam_1-to-1.nr.txt'

file1 = '/Users/biggus/Documents/Programming/WingProjects/JamesProject_REPO/MDOSX_proj/testData/test_Aa2Ag.txt'
file2 = '/Users/biggus/Documents/Programming/WingProjects/JamesProject_REPO/MDOSX_proj/testData/test_Aa2Cq.txt'
file3 = '/Users/biggus/Documents/Programming/WingProjects/JamesProject_REPO/MDOSX_proj/testData/test_Ag2Cq.txt'

oFile = '/Users/biggus/Desktop/combineOrthos.out.filtered.txt'

t1 = time()
orthoTab = combineOrthoTabs([file1,file2,file3],3)
t2 = time()

print 'merging took %s min.' % ((t2-t1)/60.0)

oFile = open(oFile, 'w')

for i in orthoTab:
    oFile.write(str(i)+'\n')