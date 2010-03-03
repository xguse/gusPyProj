

import csv
import matplotlib as mpl
#mpl.use('QtAgg')
mpl.use('TkAgg')
from matplotlib import pylab as pl

inF = '/Users/biggus/Documents/James/Data/mariangela/Normalised_expression_values_per_transcript.csv'
inF = map(lambda l: l.strip('\n').split(','), open(inF, 'rU'))


title  = None
xLab   = 'Transcripts'
yLab   = 'Fraction of Reads (Normalized)'
saveAs = '/Users/biggus/Documents/James/Data/mariangela'



def cleanData(inF,col):
    header = inF[0][col]
    tmpData = []
    for i in range(1,len(inF)):
        try:
            #datum = float(inF[i][col])
            #if datum != 0:
                #tmpData.append(datum)
            tmpData.append(float(inF[i][col]))
        except ValueError:
            print "ValueError: %s" % (inF[i][col])
    return (header,tmpData)

title,data = cleanData(inF,5)


data.sort()
data.reverse()

fig = pl.figure()
ax  = fig.add_subplot(111)

ax.plot(data)

#ax.hist(data, bins=50, range=None, normed=False, cumulative=False,
        #bottom=None, histtype='bar', align='mid',
        #orientation='vertical', rwidth=None, log=1)


ax.set_xlabel(xLab)
ax.set_ylabel(yLab)
ax.set_yscale('log')
ax.set_title(title)
pl.savefig('%s/%s_fracVsTx.pdf' % (saveAs,title))
pl.show()





        