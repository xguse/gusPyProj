from MDOSX_defs.defs import spawnOrthoGroups


nWayOrthoList = [['AAEL1','AGAP1','CPJI1'],
                 ['AAEL3','AGAP3','CPJI3'],
                 ['AAEL5','AGAP5','CPJI5'],]


o1 = '/Users/biggus/Documents/Programming/WingProjects/JamesProject_REPO/MDOSX_proj/testData/test_AgPromoters.fas'
o2 = '/Users/biggus/Documents/Programming/WingProjects/JamesProject_REPO/MDOSX_proj/testData/test_AaPromoters.fas'
o3 = '/Users/biggus/Documents/Programming/WingProjects/JamesProject_REPO/MDOSX_proj/testData/test_CqPromoters.fas'



oGroups = spawnOrthoGroups([o1,o2,o3],nWayOrthoList)
None


