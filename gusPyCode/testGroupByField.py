from JamesDefs import groupByField
import string

clusterDefinitionList = map(string.strip, open('/Users/biggus/Documents/MBGB/Rotations/James/Data/ClusterHyperGeo/testClusterDefs.txt', 'r'))


clusterDefinitionList = groupByField(clusterDefinitionList, 0)


lookInSide = clusterDefinitionList[0][0][1]

print 'yay'