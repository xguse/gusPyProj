expnFile  = map(lambda l: l.strip(), open('','rU').readlines())
namesFile = map(lambda l: l.strip(), open('','rU').readlines())
outFile   = ''

foundVectors   = []
missingVectors = []

for gene in namesFile:
    tempVector = []
    
    for vector in expnFile:
        if vector.startswith(gene):
            tempVector.append(vector)
    
    if tempVector:
        foundVectors.extend(tempVector)
    else:
        missingVectors.append(gene)

print 'The following %s genes were not fount in the vector file:\n%s' % \
      (len(missingVectors), '\n'.join(missingVectors))
        
outFile = open(outFile, 'w')

for i in foundVectors:
    outFile.write(i)

            
    

