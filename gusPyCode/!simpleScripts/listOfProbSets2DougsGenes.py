goodProbes = map(lambda l: l.strip(), open('/Users/biggus/Documents/James/Data/DoubleBloodMeal/edDBMclusters/edDBMclus4.probSets.txt', 'rU').readlines())
probe2genList = map(lambda l: l.strip('\n').split('\t'), open('/Users/biggus/Documents/James/Data/MicroArray/OMsAgExperiments/probSet2geneList.txt','rU').readlines())

prob2geneDict = {}
for i in range(0,len(probe2genList)):
    prob2geneDict[probe2genList[i][0]] = probe2genList[i][1] 

genes    = []
notFound = []

for i in range(0,len(goodProbes)):
    if goodProbes[i] in prob2geneDict:
        if prob2geneDict[goodProbes[i]] in genes:
            continue
        else:
            genes.append(prob2geneDict[goodProbes[i]])
    else:
        notFound.append(goodProbes[i])
    
     
print len(genes)
for e in genes:
    print e

print 'There were %s probes not found in ID table.' % (len(notFound))   
for n in notFound:
    print n

None
