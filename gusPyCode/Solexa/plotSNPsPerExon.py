

import csv
import matplotlib as mpl
#mpl.use('QtAgg')
mpl.use('TkAgg')
from matplotlib import pylab as pl

inF = '/Users/biggus/Documents/James/Data/Solexa/snps/RBandS.minus_rRNA.varscan.snps.combined.readcounts.filtered.snpsPerBp.txt'
inF = map(lambda l: l.strip('\n').split('\t'), open(inF, 'rU'))


title  = 'RexD B+S cov=50'
xLab   = 'Exons'
yLab   = 'SNPs/bp'
saveAs = '/Users/biggus/Documents/James/Data/Solexa/snps/'



def parsefFile(inF,col):
    tmpData = []
    for i in range(1,len(inF)):
        try:
            #datum = float(inF[i][col])
            #if datum != 0:
                #tmpData.append(datum)
            tmpData.append(float(inF[i][col]))
        except ValueError:
            print "ValueError: %s" % (inF[i][col])
    return tmpData

data = parsefFile(inF,5)



data.sort()
data.reverse()

fig = pl.figure()
ax  = fig.add_subplot(111)

ax.plot(data,linewidth=3)

#ax.hist(data, bins=50, range=None, normed=False, cumulative=False,
        #bottom=None, histtype='bar', align='mid',
        #orientation='vertical', rwidth=None, log=1)


ax.set_xlabel(xLab)
ax.set_ylabel(yLab)
#ax.set_yscale('log')
ax.set_title(title)
pl.savefig('%s/%s.pdf' % (saveAs,title))
pl.show()





        