sList = [[[1,2],[2,3],[3,4],],
         [[4,5],[5,6],[6,7],],
         [[7,8],[8,9],[9,10]]]

for i in range(len(sList)):
    for j in range(len(sList[i])):
        sList[i][j] = frozenset(sList[i][j])
    sList[i] = frozenset(sList[i])

sSet = set(sList)

tSet = set([9,10])

if tSet in sSet:
    print 'yes'
None