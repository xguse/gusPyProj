import re
from time import time

#========================= User Defined Variables =========================
inFile      = '/Users/biggus/Documents/James/Data/2KB/2kb_Sequence/2kb_Combo/2Kb_AllMosquitoes/MosqMotifs/outPut/outPut/2KBup_collapsed_upstream_ex-conserved_mosquito-motifs_nr.te.hgp.txt'

outFile = open('/Users/biggus/Documents/James/Data/2KB/2kb_Sequence/2kb_Combo/2Kb_AllMosquitoes/MosqMotifs/outPut/outPut/2KBup_collapsed_upstream_ex-conserved_mosquito-motifs_nr.te.hgp.ngt0.txt','w')

fieldIndex = 2
#===============================================

inFile = map(lambda line: line.strip(), open(inFile, 'rU').readlines())
len_inFile = len(inFile)

culledLines =[]

regEx = re.compile('protein_coding', re.IGNORECASE)
for item in inFile:
    
    itemFields = item.split('\t')
    if itemFields[4] != '0':    #  regEx.search(itemFields[fieldIndex])
        
        culledLines.append(item+'\n')
        


len_culledLines = len(culledLines)

print '%s Lines were retrieved\n%s lines were present originaly\n%i lines were removed' % (str(len_culledLines),str(len_inFile), len_inFile - len_culledLines)

outFile.writelines(culledLines)