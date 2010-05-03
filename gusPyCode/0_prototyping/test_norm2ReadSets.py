from gusPyCode.defs.prototypes import norm2ReadSets

l1 = [1,2,3,4,5,6,7,8,9]
l2 = [2,5,5,100,60,7,6,9,9]


nrmd1,nrmd2 = norm2ReadSets(l2,l1)

print nrmd1
print nrmd2