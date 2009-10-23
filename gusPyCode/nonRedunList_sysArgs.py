import sys
import re



#print 'sys args = %s' % (len(sys.argv))
assert len(sys.argv) == 2, 'Usage: nonRedunList_sysArgs.py listFile'

############### list of user variables ####################

# file to import from
rFile = open(sys.argv[1],'r').readlines()


# file to write output to
xtn = re.compile('\.\w{1,20}?.$')
oldXtn = xtn.findall(sys.argv[1])
if oldXtn:
    outFile = xtn.sub('.nr%s' % (oldXtn[0]),sys.argv[1])
else:
    outFile = sys.argv[1]+'.nr'

# open and read data into var



# invoke list to hold results
nrList = []
repeatList = []

# read one line from file at a time and:
# -strip the trailing \n
# -ask if it is in the results list
# -if no put it there, if yes go to next line

for line in rFile:
    line = line.rstrip()
    if line not in nrList:
        nrList.append(line)
    else:
        if line not in repeatList:
            repeatList.append(line)
        


#open outFile and write data/close file
nrFile = open(outFile, 'w')
nrFile.write('\n'.join(nrList))

print str(len(nrList))+' non-repetative sequences written. '+str(len(rFile)-len(nrList))+' items excluded.'
print 'The following is a non-redundant list of those records that were repeated in original file:'
for everyOne in repeatList:
    print everyOne