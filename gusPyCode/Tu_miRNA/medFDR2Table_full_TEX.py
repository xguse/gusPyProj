"""
Read in medFDR outFile and output .tex version of table for paper with only the combined data.
"""

import sys
from gusPyCode.defs.JamesDefs import groupByField_silent
from gusPyCode.defs.JamesDefs import DotDict
if len(sys.argv) != 3:
    exit('USAGE: python %s inFile outFile' % (sys.argv[0].split('/')[-1]))

    
# ----- DEFS -----
def getAGAPList(strOfTupReps):
    rSet = set()
    tups = eval(strOfTupReps)
    
    for tup in tups:
        for i in tup:
            if i.startswith('AGAP'):
                rSet.add(i)
    return sorted(list(rSet))
    
    
    
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
            data.numTot[orthoType] = len(getAGAPList(line[5]))
            data.FDRTot[orthoType] = line[3]
            for i in miR:
                if i[2].startswith('orthoType_%s' % (orthoType)):
                    if i[1]+'_fdr' in data:
                        data[i[1]+'_fdr'][orthoType] = i[4]
            
    outTmp = []
    oClassMap = {1:'I',
                 2:'II',
                 3:'III'}
    for i in range(len(data.orthoTypes)):
        data.orthoTypes.sort()
        if i == 0:
            
            outTmp.append('\multirow{%s}{*}{%s} & %s & %s & %.3f & %s & %s & %s & %s & %s\\\\\n'%\
                          (str(len(data.orthoTypes)),
                           data.name,
                           oClassMap[data.orthoTypes[i]],
                           data.numTot[data.orthoTypes[i]],
                           float(data.FDRTot[data.orthoTypes[i]]),
                           data.A1_to_m7_fdr[data.orthoTypes[i]],
                           data.A1_to_m8_fdr[data.orthoTypes[i]],
                           data.m2_to_m7_fdr[data.orthoTypes[i]],
                           data.m2_to_m8_fdr[data.orthoTypes[i]],
                           data.m3_to_m8_fdr[data.orthoTypes[i]]))
        else:
            outTmp.append('& %s & %s & %.3f & %s & %s & %s & %s & %s\\\\\n'% \
                          (oClassMap[data.orthoTypes[i]],
                           data.numTot[data.orthoTypes[i]],
                           float(data.FDRTot[data.orthoTypes[i]]),
                           data.A1_to_m7_fdr[data.orthoTypes[i]],
                           data.A1_to_m8_fdr[data.orthoTypes[i]],
                           data.m2_to_m7_fdr[data.orthoTypes[i]],
                           data.m2_to_m8_fdr[data.orthoTypes[i]],
                           data.m3_to_m8_fdr[data.orthoTypes[i]]))
        

    outTmp.append('\\hline\n')
    outLines.extend(outTmp)
    
    
outFile = open(sys.argv[2], 'w')
outFile.write('\\begin{center}\n\\begin{tabular}{lcccccccc}')
outFile.write('\\hline\\hline\n')
outFile.write('miRNA & Success Class & Combined Targets & Combined FDR & A1-m7 & A1-m8 & m2-m7 & m2-m8 & m3-m8 \\\\ [0.5ex] \\hline\n')
for line in outLines:
    outFile.write(line)
outFile.write('\\end{tabular}')
    
print "I'm Done!\nSee file: %s" % (sys.argv[2])
        