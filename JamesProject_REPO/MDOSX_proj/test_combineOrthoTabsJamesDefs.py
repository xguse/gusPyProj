from JamesDefs import combineOrthoTabs
from pprint import pprint
from time import time


files = ['/Users/biggus/Documents/James/Data/OrthologDefs/nrOrthos/!nrOrthoDefs/Aaeg_Agam_1-to-1.nr.txt',
         '/Users/biggus/Documents/James/Data/OrthologDefs/nrOrthos/!nrOrthoDefs/Aedes_Culex_1-to-1.nr.txt',
         '/Users/biggus/Documents/James/Data/OrthologDefs/nrOrthos/!nrOrthoDefs/Culex_Agam_1-to-1.nr.txt',]



#['/Users/biggus/Documents/Programming/WingProjects/JamesProject_REPO/MDOSX_proj/testData/test_Aa2Ag.txt',
         #'/Users/biggus/Documents/Programming/WingProjects/JamesProject_REPO/MDOSX_proj/testData/test_Aa2Cq.txt',
         #'/Users/biggus/Documents/Programming/WingProjects/JamesProject_REPO/MDOSX_proj/testData/test_Ag2Cq.txt',
         #'/Users/biggus/Documents/Programming/WingProjects/JamesProject_REPO/MDOSX_proj/testData/test_Ag2Dm.txt',
         #'/Users/biggus/Documents/Programming/WingProjects/JamesProject_REPO/MDOSX_proj/testData/test_Aa2Dm.txt',
         #'/Users/biggus/Documents/Programming/WingProjects/JamesProject_REPO/MDOSX_proj/testData/test_Cq2Dm.txt',]
# list-ify file data
for i in range(len(files)):
    files[i] = map(lambda l: l.strip().split(), open(files[i], 'rU').readlines())


oFile = '/Users/biggus/Desktop/combineOrthosJamesDefs.out.filtered.txt'

t1 = time()
orthoTab = combineOrthoTabs(files)
t2 = time()

print 'merging took %s min.' % ((t2-t1)/60.0)


oFile = open(oFile, 'w')

for i in orthoTab:
    oFile.write(str(i)+'\n')