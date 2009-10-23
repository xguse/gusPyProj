from MDOSX_defs.defs import buildCoRegdOrthoGroups


nWayOrthoList = [['A1','B3','C5'],
                 ['A3','B7','C1'],
                 ['A8','B1','C4'],]


coRegdGeneListFile = '/Users/biggus/Documents/Programming/WingProjects/JamesProject_REPO/MDOSX_proj/testData/test_coRegdGeneList.txt'

coRegdOrthoGroups = buildCoRegdOrthoGroups(nWayOrthoList, coRegdGeneListFile)
None


