from matplotlib import pyplot as plt
from matplotlib import pylab as plb
from gusPyCode.defs import JamesDefs
import miRNA_targeting as miTrgt
import sys
from pprint import pprint

#orthoDefs = [['Ag_oA','Ae_oA','Cq_oA'],
             #['Ag_oB','Ae_oB','Cq_oB'],
             #['Ag_oC','Ae_oC','Cq_oC'],
             #['Ag_oD','Ae_oD','Cq_oD']]

#miRNAs = {'aga-bantam':'UGAGAUCACUUUGAAAGCUGAUU',
          #'aga-let-7':'UGAGGUAGUUGGUUGUAUAGU',
          #'aga-miR-1':'UGGAAUGUAAAGAAGUAUGGAG',
          #'aga-miR-10':'ACCCUGUAGAUCCGAAUUUGU',
          #'aga-miR-100':'AACCCGUAGAUCCGAACUUGUG'}
print '-- -- -- -- --\n\nLoading files...'
orthoPath  = '/Users/biggus/Documents/James/Data/OrthologDefs/nrOrthos/!nrOrthoDefs/3wayOrthos.combineOrthosJamesDefs.out.filtered.noInfers.txt'
seqPaths   = ['/Users/biggus/Documents/James/Data/Tu_miRNA/Fastas/Aa_500afterCoding.usuable.stpCdn.fas',
              '/Users/biggus/Documents/James/Data/Tu_miRNA/Fastas/Ag_500afterCoding.usuable.stpCdn.fas',
              '/Users/biggus/Documents/James/Data/Tu_miRNA/Fastas/Cq_500afterCoding.newCoords.usuable.stpCdn.fas']          
miRNA_Path = '/Users/biggus/Documents/James/Data/Tu_miRNA/miRNAs/miRBase/mature.aga.fa'

          
seqDict        = miTrgt.loadSeqs(seqPaths)
orthoRelations = miTrgt.loadOrthos(orthoPath, seqDict)
miRNAs         = miTrgt.loadMiRNAs(miRNA_Path)

print 'Picking random data...'
randOrthos = orthoRelations #JamesDefs.randFromList_noRplcMulti(orthoRelations,1000)
randMiRNAs = ['aga-miR-11'] #JamesDefs.randFromList_noRplcMulti(miRNAs.keys(),1)
orthoSeqs  = miTrgt.filterOrthoSeqs(seqDict,randOrthos)


seenSeeds = set()
miR_matches = {}




print 'Initializing matchVersions...'

for m in randMiRNAs:
    seed = miTrgt.miRNA(miRNAs[m],seenSeeds,orthoRelations=randOrthos,name=m)
    miR_matches[seed.name]= seed

print 'Initializing ctrls...'
for m in miR_matches:
    miR_matches[m].buildCtrlsFromProSeed(seenSeeds,5)


for m in miR_matches:
    print miR_matches[m].name
    for sVer in miR_matches[m].matchVersions:
        print '%s: %s' % (sVer, miR_matches[m].matchVersions[sVer])
    print '- '*5
print 'Tallying hits...'
for m in miR_matches:
    miR_matches[m].tallyHits(orthoSeqs) 
print 'Counting hits in orthos...'
for m in miR_matches:
    miR_matches[m].countHitsInOrthos() 



print 'Drawing BoxPlots...'


# Figure options
xlabel = 'matchType'
ylabel = 'counts'
title ='%s in all orthoRelations' % (m)
saveToFile = 'countData.%s.boxPlot.pdf' % (m)
dpi = 1200

seedTypes = sorted(miTrgt._seedModels)

boxData = []
labelTexts = []
for m in sorted(miR_matches.keys()):
    for seedType in seedTypes:
        for i in range(1,len(miR_matches[m].ctrlCounts[seedType][0])):
            boxData.append([x[i] for x in miR_matches[m].ctrlCounts[seedType]])
            labelTexts.append('%s:%s:%s' % (m,seedType,i))

realData = []
rLabels  = []
for m in sorted(miR_matches.keys()):
    for seedType in seedTypes:
        realData.extend(miR_matches[m].matchCounts[seedType][1:])
        for i in range(1,4):
            rLabels.append('%s:%s:%s' % (m,seedType,i))

#pprint(realData)       
           
    
pos  = range(1,len(boxData)+1)
#posR = range(pos[-1]+1,pos[])

fig = plb.figure()
fig.suptitle(title)
ax = fig.add_subplot(1,1,1)

bpInfo = ax.boxplot(boxData, notch=0, sym='x', vert=1, whis=1.5,positions=None, widths=None,)

ax2 = fig.add_subplot(1,1,1)
barInfo = ax2.bar(pos,realData,color='g',align='center')

    
ax.set_ylabel(ylabel)
ax2.set_xlabel(xlabel)

ax2.set_xticks(pos)
ax2.set_xticklabels(labelTexts)
fig.autofmt_xdate()
fig.set_size_inches(50,20)

#plt.axis([0, pos[-1], 0, max(max(max(boxData)))])
plt.grid(False)
#plt.legend(loc=4)
if saveToFile:
    plt.savefig(saveToFile,dpi=dpi)
plt.show()

    
None
