"""
Read in medFDR outFile and output tsv version of table for paper.
"""

import sys
from gusPyCode.defs.JamesDefs import groupByField_silent
from gusPyCode.defs.JamesDefs import DotDict
if len(sys.argv) != 3:
    exit('USAGE: python %s inFile outFile' % (sys.argv[0].split('/')[-1]))

    
lines = open(sys.argv[1],'rU').readlines()
outLines = []

groupedLines = groupByField_silent(lines,0,sep=' : ')
groupedLines.sort(key=lambda x: x[0][0])
for miR in groupedLines:
    data = DotDict({'name':miR[0][0],
                       'orthoTypes':[],
                       'numTot':[None]*4,
                       'FDRTot':[None]*4,
                       'm3_to_m8_fdr':[None]*4,
                       'm2_to_m7_fdr':[None]*4,
                       'A1_to_m7_fdr':[None]*4,
                       'm2_to_m8_fdr':[None]*4,
                       'A1_to_m8_fdr':[None]*4})
    
    for line in miR:
        if line[1].startswith("allPassedSeedsFor_"):
            orthoType = int(line[1][-1])
            data.orthoTypes.append(orthoType)
            data.numTot[orthoType] = line[2]
            data.FDRTot[orthoType] = line[3]
            for i in miR:
                if i[2].startswith('orthoType_%s' % (orthoType)):
                    if i[1]+'_fdr' in data:
                        data[i[1]+'_fdr'][orthoType] = i[4]
            
    outTmp = []        
    for i in range(len(data.orthoTypes)):
        data.orthoTypes.sort()
        outTmp.append('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % \
                      (data.name,
                       data.orthoTypes[i],
                       data.numTot[data.orthoTypes[i]],
                       data.FDRTot[data.orthoTypes[i]],
                       data.A1_to_m7_fdr[data.orthoTypes[i]],
                       data.A1_to_m8_fdr[data.orthoTypes[i]],
                       data.m2_to_m7_fdr[data.orthoTypes[i]],
                       data.m2_to_m8_fdr[data.orthoTypes[i]],
                       data.m3_to_m8_fdr[data.orthoTypes[i]]))
    outLines.extend(outTmp)
    
    
outFile = open(sys.argv[2], 'w')
outFile.write('miRNA\t--\tCombined Targets\tCombined FDR\tA1-m7\tA1-m8\tm2-m7\tm2-m8\tm3-m8\n')
for line in outLines:
    outFile.write(line)
    
print "I'm Done!\nSee file: %s" % (sys.argv[2])
        