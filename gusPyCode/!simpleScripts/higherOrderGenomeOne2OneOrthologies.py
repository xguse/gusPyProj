from gusPyCode.defs.xpermutations import *



def forNgenomes(numOfGenomes):
    x = 2
    results = []
    while x <= numOfGenomes:
        temp = []
        for uc in xuniqueCombinations(range(numOfGenomes),x):
            temp.append(uc)
            
        results.append(temp[:])
        x += 1 
    return results


n0 = 3
n  = 6
answers = []

for numOfGenomes in range(n0,n+1):
    answers.append(forNgenomes(numOfGenomes))
    

for i in answers:
    print "For %s genomes:" % (len(i)+1)
    c = 2
    for j in i:
        print "Combos of %s:\t%s" % (c,len(j),) #str(j))
        #print "%s\t%s" % (c,len(j),)  #  Just the numbers
        c+=1

None
