from defs_microArray import geneVectors


l = ['one\t2.4\t4.3\t8.3','two\t2.4\t4.3\t8.3','three\t2.4\t4.3\t8.3']

vectors = geneVectors(l)


for i in vectors.vectors:
    print '%s:%s' % (i,vectors.vectors[i])
    

None
