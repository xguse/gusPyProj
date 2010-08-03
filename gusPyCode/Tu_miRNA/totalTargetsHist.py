"""Hists for total targets per miR.
"""

#import matplotlib as mpl
#mpl.use('WXAgg')
from matplotlib import pylab as pl
from gusPyCode.defs.mpl_custom import setTickSizes



path = '/Users/biggus/Documents/James/Data/Tu_miRNA/ResultsMiRNA_Targeting/2009_11_11/2009_11_11.AGAP.stats.100.medFDRmeth.cls3TrgtHist.txt'

data = map(lambda l: int(l.strip()), open(path, 'rU'))

saveAs = '/Users/biggus/Documents/James/Data/Tu_miRNA/ResultsMiRNA_Targeting/2009_11_11/2009_11_11.AGAP.stats.100.medFDRmeth.cls3TrgtHist.ntit.pdf'
title   = ''#'Class III'
xAxis   = 'Predicted Targets'
yAxis1  = 'miRNAs'
#yAxis2  = 'Proportion of Total'

fig = pl.figure()
ax1 = fig.add_subplot(111)
pl.ylim(0,40)
#ax2 = ax1.twinx()
ax1.set_title(title)
ax1.set_xlabel(xAxis,fontsize='xx-large')
ax1.set_ylabel(yAxis1,fontsize='xx-large')
#ax2.set_ylabel(yAxis2)

n,b,p = ax1.hist(data,bins=20, range=None, normed=0,
                 weights=None, cumulative=1, bottom=None,
                 histtype='stepfilled', align='mid', orientation='vertical',
                 rwidth=None, log=False,alpha=1)

#eightyPcnt = round(n[-1]*0.8)
#pl.axhspan(eightyPcnt, eightyPcnt, xmin=0, xmax=1, color='red')

#ax2.hist(data, bins=50, range=None, normed=1,
         #weights=None, cumulative=1, bottom=None,
         #histtype='stepfilled', align='mid', orientation='vertical',
         #rwidth=None, log=False, alpha=0.5, facecolor='blue')

pl.xlim(0,max(data))
setTickSizes(ax1,18)

pl.savefig(saveAs)
pl.show()
