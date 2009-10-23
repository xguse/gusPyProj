from random import choice
import cPickle
import miRNA_targeting as miRT
from matplotlib import pylab


print '\n\n'

outPath = '/Users/biggus/Documents/James/Data/Tu_miRNA/ResultsMiRNA_Targeting/seedCtrlHistBarPlots_2009_08_31/'
pklPath = '/Users/biggus/Documents/James/Data/Tu_miRNA/ResultsMiRNA_Targeting/2009_08_31.seedMatches.100ctrls.pkl'
data    = cPickle.load(open(pklPath,'rU'))

# calc the various centers and spreads
items = ['aga-miR-1890',
         'aga-miR-9c',
         'aga-miR-13b',
         'aga-miR-133',
         'aga-miR-11',
         'aga-miR-988',
         'aga-miR-278',
         'aga-miR-965',]

miR_matches = {}
for i in items:
    miR_matches[i] = data['miR_matches'][i]

#figs = []
    
# Plot hist along x-axis and a horizontal boxPlot
def drawEm(miRNA):
    for seedType in sorted(miRNA.matchVersions):
        for orthoType in range(1,4):
            f = pylab.figure()
            tit = '%s %s %s' % (miRNA.name,seedType,orthoType)
            f.suptitle(tit)
            a = f.add_subplot(111)
            data = [x[orthoType] for x in miRNA.ctrlCounts[seedType]]
            hData = a.hist(data, bins=20,facecolor='g', alpha=0.6)
            # set x and y limits based on the hist data
            xMin,xMax = int(min(hData[1])-1),int(max(hData[1])+1)
            yMin,yMax = 0,int(max(hData[0])*1.05)
            a.set_xlim(xMin,xMax) 
            a.boxplot(data,notch=1,widths=yMax/4, vert=0,positions=([yMax/2]))
            a.set_ylim(yMin,yMax)
            a.set_yticks(range(yMin,yMax,yMax//5))
            pylab.savefig('%s%s.pdf'%(outPath,tit.replace(' ','_')))
            #figs.append(f)

#miRNA = miR_matches[items[6]]

#drawEm(miRNA)

for miRNA in items:
    print 'Drawing %s' % (miRNA)
    drawEm(miR_matches[miRNA])

#pylab.show()

print 'Im Done!'
    

    
    




