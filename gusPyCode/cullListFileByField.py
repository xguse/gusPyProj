import re
from time import time

#========================= User Defined Variables =========================
inFile      = '/Users/biggus/Documents/James/Data/Tu_miRNA/Cq_CpipJ1-2_GeneTranscrExon_112408.txt'

outFile = open('/Users/biggus/Documents/James/Data/Tu_miRNA/Cq_CpipJ1-2_GeneTranscrExon_112408.noBlanks.txt','w')

fieldIndex = 2

regEx = ''

#===============================================

inFile = map(lambda line: line.strip(), open(inFile, 'rU').readlines())
len_inFile = len(inFile)

culledLines =[]

regEx = re.compile(regEx, re.IGNORECASE)
for item in inFile:
    
    itemFields = item.split('\t')
    if itemFields[5] and itemFields[6] != '':    #  regEx.search(itemFields[fieldIndex])
        
        culledLines.append(item+'\n')
        


len_culledLines = len(culledLines)

print '%s Lines were retrieved\n%s lines were present originaly\n%i lines were removed' % (str(len_culledLines),str(len_inFile), len_inFile - len_culledLines)

outFile.writelines(culledLines)