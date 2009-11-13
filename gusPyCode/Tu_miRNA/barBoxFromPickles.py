import cPickle
from matplotlib import pylab as plb
from gusPyCode.defs import miRNA_targeting as miRT

print '\n\n'

outFile = '/Users/biggus/Documents/James/Data/Tu_miRNA/ResultsMiRNA_Targeting/2009_11_11/2009_11_11.AGAP.seedMatches.300mvCtrls.stats.txt'
outFile = open(outFile, 'w')

pklPath = '/Users/biggus/Documents/James/Data/Tu_miRNA/ResultsMiRNA_Targeting/2009_11_11/2009_11_11.AGAP.seedMatches.300mvCtrls.storeEvents.pkl'
data    = cPickle.load(open(pklPath,'rU'))

dirPath = '/Users/biggus/Documents/James/Data/Tu_miRNA/ResultsMiRNA_Targeting/2009_11_11/barBoxes/'
genome  = 'AGAP'

# calc the various centers and spreads
items = ['aga-miR-1890',
         'aga-miR-9c',
         'aga-miR-13b',
         'aga-miR-133',
         'aga-miR-11',
         'aga-miR-988',
         'aga-miR-278',
         'aga-miR-965',]

#miR_matches = {}
#for i in items:
    #miR_matches[i] = data[i]

miR_matches = data
    
print 'Drawing BoxPlots...'




seedTypes = sorted(miRT._seedModels)

for m in miR_matches: 
    
    # Figure options
    xlabel = 'matchType'
    ylabel = 'counts'
    title ='%s (300 mvCtrls)' % (m)
    saveToFile = dirPath+'2009_11_11.AGAP.%s.boxPlot.300mvCtrls.pdf' % (m)
    dpi = 1200
    
    boxData = []
    labelTexts = []
    for seedType in seedTypes:
        for i in range(1,len(miR_matches[m].ctrlEvents[seedType][0])):
            boxData.append([len(miRT.filterToken(x[i],'AGAP')) for x in miR_matches[m].ctrlEvents[seedType]])
            labelTexts.append('%s : %s' % (seedType,i))

    realData = []
    rLabels  = []
    for seedType in seedTypes:
        realData.extend([len(miRT.filterToken(x,'AGAP')) for x in miR_matches[m].matchEvents[seedType][2:]])
        for i in range(1,4):
            rLabels.append('%s : %s' % (seedType,i))
    
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
    
    
    #plt.axis([0, pos[-1], 0, max(max(max(boxData)))])
    plb.grid(False)
    #plt.legend(loc=4)
    if saveToFile:
        plb.savefig(saveToFile,dpi=dpi)
#plb.show()


print 'Im Done.'