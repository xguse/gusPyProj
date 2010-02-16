"""
Read in medFDR outFile and output CoOccurring miRs.
"""

import sys
from gusPyCode.defs.JamesDefs import groupByField_silent
from gusPyCode.defs.JamesDefs import DotDict
from gusPyCode.defs.JamesDefs import Bag
from gusPyCode.defs.xpermutations import xuniqueCombinations

if len(sys.argv) != 3:
    exit('USAGE: python %s inFile outFile' % (sys.argv[0].split('/')[-1]))
# ----- DEFS -----
def getAGAP(strOfTupReps):
    rSet = set()
    tups = eval(strOfTupReps)
    
    for tup in tups:
        for i in tup:
            if i.startswith('AGAP'):
                rSet.add(i)
    return rSet
            
    
    
lines = open(sys.argv[1],'rU').readlines()
outLines = []

groupedLines = groupByField_silent(lines,0,sep=' : ')
groupedLines.sort(key=lambda x: x[0][0])

mirDict = Bag()

for miR in groupedLines:
    data = DotDict({'name':miR[0][0],
                    'orthoTypes':[],
                    'numTot':[None]*4,
                    'AGAPgenes':[None]*4})
    
    for line in miR:
        if line[1].startswith("allPassedSeedsFor_"):
            orthoType = int(line[1][-1])
            data.orthoTypes.append(orthoType)
            data.numTot[orthoType] = len(getAGAP(line[5]))
            data.AGAPgenes[orthoType] = getAGAP(line[5])
            
            #for i in miR:
                #if i[2].startswith('orthoType_%s' % (orthoType)):
                    #if i[1]+'_fdr' in data:
                        #data[i[1]+'_fdr'][orthoType] = i[4]
    
    mirDict[data.name] = data

# ---- Work out combinations ----
mirCombos = sorted([sorted(x) for x in xuniqueCombinations(mirDict.keys(),2)])
l = len(mirCombos)
outTmp = []        
for mCombo in mirCombos:
    for i in range(2,4):
        if (mirDict[mCombo[0]].numTot[i] == None) or (mirDict[mCombo[1]].numTot[i] == None):
            pass
        else:
            cmbo  = ':'.join(mCombo)
            clas  = 'Class %s' % (i)
            eatot = '%s:%s' % (mirDict[mCombo[0]].numTot[i],mirDict[mCombo[1]].numTot[i])
            c1Set = mirDict[mCombo[0]].AGAPgenes[i]
            c2Set = mirDict[mCombo[1]].AGAPgenes[i]
            inTot = '%s' % (len(c1Set.intersection(c2Set)))
            inter = '%s' % (sorted(list(c1Set.intersection(c2Set))))
            
            outTmp.append('%s\t%s\t%s\tnum_%s\t%s\n' % \
                          (cmbo,
                           clas,
                           eatot,
                           inTot,
                           inter))
outLines.extend(outTmp)
    
    
outFile = open(sys.argv[2], 'w')
outFile.write('miRNA Combo\tClass\tTotal Targets\tNumber Co-occurring\tGenes\n')
for line in outLines:
    outFile.write(line)
    
print "I'm Done!\nSee file: %s" % (sys.argv[2])
        