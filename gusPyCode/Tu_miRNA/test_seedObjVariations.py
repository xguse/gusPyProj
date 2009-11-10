from matplotlib import pylab as plb
from gusPyCode.defs import JamesDefs
import miRNA_targeting as miTrgt
import sys



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

print 'Filtering orthos...'
orthoSeqs  = miTrgt.filterOrthoSeqs(seqDict,orthoRelations)


seenSeeds = set()
miR_matches = {}

saveObj = {'orthoRelations':orthoRelations,
           'seenSeeds':seenSeeds,
           'miR_matches':miR_matches}


print 'Initializing matchVersions...'

for m in miRNAs:
    seed = miTrgt.seedMatches(miRNAs[m],seenSeeds,orthoRelations=orthoRelations,name=m)
    miR_matches[seed.name]= seed

print 'Initializing ctrls...'
# choose one rand miRNA to make ctrls
randMiRNA = JamesDefs.randFromList_noRplcMulti(miRNAs.keys(),1)[0]
miR_matches[randMiRNA].buildCtrlsFromMatchVers(seenSeeds,30)
randMiRNA = miR_matches[randMiRNA]


#for m in miR_matches:
    #print m.name
    #for sVer in m.matchVersions:
        #print '%s: %s' % (sVer, m.matchVersions[sVer])
    #print '- '*5
print 'Tallying hits...'
randMiRNA.tallyHits(orthoSeqs) 
print 'Counting hits in orthos...'
randMiRNA.countHitsInOrthos() 


exit('No Plotting today!')
print 'Drawing BoxPlots...'


# Figure options
xlabel = 'matchType'
ylabel = 'counts'
title ='%s in all orthoRelations' % (randMiRNA.name)
saveToFile = '/Users/biggus/Documents/James/Data/Tu_miRNA/ResultsMiRNA_Targeting/countData.%s.boxPlot.pdf' % (randMiRNA.name)
dpi = 1200

seedTypes = sorted(miTrgt._seedModels)

boxData = []
labelTexts = []

for seedType in seedTypes:
    for i in range(1,len(randMiRNA.ctrlCounts[seedType][0])):
        boxData.append([x[i] for x in randMiRNA.ctrlCounts[seedType]])
        labelTexts.append('%s:%s:%s' % (randMiRNA.name,seedType,i))

realData = []
rLabels  = []
for seedType in seedTypes:
    realData.extend(randMiRNA.matchCounts[seedType][1:])
    for i in range(1,4):
        rLabels.append('%s:%s:%s' % (randMiRNA.name,seedType,i))

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
#fig.set_size_inches(50,20)

#plt.axis([0, pos[-1], 0, max(max(max(boxData)))])
plt.grid(False)
#plt.legend(loc=4)
if saveToFile:
    plt.savefig(saveToFile,dpi=dpi)
plt.show()

    
None
