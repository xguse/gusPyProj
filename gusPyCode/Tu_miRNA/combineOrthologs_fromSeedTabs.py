import sys
import JamesDefs


#print 'sys args = %s' % (len(sys.argv))
assert len(sys.argv) == 5, 'Usage: combineOrthologs_fromSeedTabs.py tabFile1 tabfile2 orthoDefs outfile'

tabFile1  = map(lambda line: line.strip(), open(sys.argv[1],'rU').readlines())
tabFile2  = map(lambda line: line.strip(), open(sys.argv[2],'rU').readlines())
orthoDefs = map(lambda line: line.strip(), open(sys.argv[3],'rU').readlines())
oFile     = sys.argv[4]


# affirm that column titles match in each tabFile and save the info
assert tabFile1[0] == tabFile2[0], 'Column titles do not match between source files.'
columnTitles = tabFile1[0]

# cleanse commented lines from both lists
tabFile1 = JamesDefs.removeCommentLines(tabFile1,'#')
tabFile2 = JamesDefs.removeCommentLines(tabFile2,'#')


    
# create one dict from tabFile1&2
combinedDict = {}
for line in tabFile1:
    fields = line.split('\t',1)
    combinedDict[fields[0]]=fields[1]

for line in tabFile2:
    fields = line.split('\t',1)
    combinedDict[fields[0]]=fields[1]
    
# write the new list
outList = ['#%s\n%s\n' % (' '.join(sys.argv),columnTitles)]
genesInDict = combinedDict.keys()
for orthoPair in orthoDefs:
    orthoPair = orthoPair.split('\t')
    if orthoPair[0] in genesInDict and orthoPair[1] in genesInDict:
        outList.append('%s\t%s\n' % (orthoPair[0],combinedDict[orthoPair[0]]))
        outList.append('%s\t%s\n' % (orthoPair[1],combinedDict[orthoPair[1]]))
        outList.append('--\n')
    else:
        print 'One or both of %s were not found in combinedDict.' % (orthoPair)


oFile = open(oFile,'w')
oFile.writelines(outList)

print "Done."
