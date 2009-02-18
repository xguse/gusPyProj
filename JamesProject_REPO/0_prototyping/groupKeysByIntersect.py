from sets import Set
from pprint import pprint
from time import time


testFile = '/Users/biggus/Documents/James/Data/2KB/2kb_Sequence/2kb_Anopheles/Modules/subSet_dictOfMotifPairsByAGAP_Anopheles.txt'

# construct Dict from file
testFile = map(lambda line: line.strip() , open(testFile, 'rU').readlines())

testDict = {}

for l in range(0, len(testFile)):
    x = testFile[l].split('\t')
    testDict[x[0]] = Set(x[1:])
    
    
# group each key with the top 3 other keys based on intersection of their sets
keys = testDict.keys()
listOfKeyLists = []

t1 = time()
for K in keys:
    matches = []
    for k in keys:
        intersect = len(testDict[K].intersection(testDict[k]))
        matches.append([intersect,k])
    matches.sort()
    matches.reverse()
    listOfKeyLists.append([K,matches[0:4]])
    
t2 = time()

print 'Comparison of %s keys took %s seconds.' % (len(keys), t2-t1)
pprint(listOfKeyLists)
                          
