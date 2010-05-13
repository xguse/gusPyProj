

#import csv
#import matplotlib as mpl
##mpl.use('QtAgg')
#mpl.use('TkAgg')

from matplotlib import pylab as pl

inF = '/Users/biggus/Documents/James/Writings_Talks/Manuscripts/AaTxDelta_RNAseq/2010_05_10/Manuscript/LBvLS_txReadCounts.txt'
inF = map(lambda l: l.strip('\n').split('\t'), open(inF, 'rU'))


title   = 'Sugarfed'
xLab    = 'Transcripts'
yLab    = 'Normalized Read Counts Mapped to Transcript'
saveAs  = '/Users/biggus/Documents/James/Writings_Talks/Manuscripts/AaTxDelta_RNAseq/2010_05_10/Manuscript'
xLim    =(0,15000)
yLim    =(1e-1,1e6)

toPlot = 1

def norm2ReadSets(list1, list2):
    """Takes two lists -> Returns two Normalized Lists.
    Normalizes Read counts in each index based
    on the list with the smallest total (this 
    list is not changed). Indexes in list with 
    larger total are multiplied by 
    (float(sum(smaller))/sum(larger)) to scale the
    read counts to those expected if both total reads
    equaled sum(smaller)."""
    if not len(list1) == len(list2):
        raise Exception, "ERROR: len(list1) != len(list2)."
    
    # ++ float-ify lists ++
    for i in range(len(list1)):
        list1[i] = float(list1[i])
        list2[i] = float(list2[i]) 
    
    totL1,totL2 = sum(list1),sum(list2)
    
    normd1,normd2 = [],[]
    
    if totL1 < totL2:
        for i in range(len(list1)):
            normd1.append(list1[i])
            normd2.append(list2[i]*(totL1/totL2))
    else:
        for i in range(len(list1)):
            normd2.append(list2[i])
            normd1.append(list1[i]*(totL2/totL1))
    
    return normd1,normd2



data = norm2ReadSets([x[1] for x in inF[1:]], [x[2] for x in inF[1:]])



data[toPlot].sort()
data[toPlot].reverse()

fig = pl.figure()
ax  = fig.add_subplot(111)

ax.plot(data[toPlot])

#ax.hist(data, bins=50, range=None, normed=False, cumulative=False,
        #bottom=None, histtype='bar', align='mid',
        #orientation='vertical', rwidth=None, log=1)


ax.set_xlabel(xLab,fontsize='large')
ax.set_ylabel(yLab,fontsize='large')
ax.set_yscale('log')
ax.set_title(title,fontsize='large')
ax.set_xlim(xLim)
ax.set_ylim(yLim)
ax.minorticks_on()
pl.savefig('%s/%s_normdReadsVsTx.pdf' % (saveAs,title))
pl.show()





        