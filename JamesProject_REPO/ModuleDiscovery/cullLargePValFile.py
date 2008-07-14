import re
from time import time

#========================= User Defined Variables =========================
pValsFile      = '/Users/biggus/Documents/James/Data/2KB/2kb_Sequence/2kb_Combo/2Kb_AllMosquitoes/2KBup_9sConservedInAaegAgamCulex.collapsedOnPerf.pVals.txt'

outFile = open('/Users/biggus/Documents/James/Data/2KB/2kb_Sequence/2kb_Combo/2Kb_AllMosquitoes/2KBup_9sConservedInAaegAgamCulex.collapsedOnPerf.culled.pVals.txt','w')
#===============================================

pValsFile = map(lambda line: line.strip(), open(pValsFile, 'rU').readlines())
len_pValsFile = len(pValsFile)

culledCombos =[]

regEx = re.compile('3h_up', re.IGNORECASE)

for item in pValsFile:
    
    itemFields = item.split('\t')
    
    
    if float(itemFields[5]) <= 0.05:    # OR regEx.search(itemFields[0])
        
        culledCombos.append(item+'\n')
        


len_culledCombos = len(culledCombos)

print '%s combos were retrieved from pValsFile which had a total of %s.' % (str(len_culledCombos),str(len_pValsFile))

outFile.writelines(culledCombos)