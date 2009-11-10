from gusPyCode.defs import JamesDefs
import string

delimitedList = open('/Users/biggus/Documents/MBGB/Rotations/James/Data/Sequence/Culex/Culex_Exon_Location.txt', 'r').readlines()

delimitedList = map(string.strip, delimitedList)

newList = JamesDefs.reOrderDelimitedList(delimitedList, '\t', [7,1,8,4,6,2,3,5,0])


outFile = open('/Users/biggus/Documents/MBGB/Rotations/James/Data/Sequence/Culex/Culex_Exon_Location_Reordered.txt','w')

for rec in newList:
    outFile.write(rec+'\n')


print 'Yay'