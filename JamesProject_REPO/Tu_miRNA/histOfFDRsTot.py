from matplotlib import pylab as plb

"""
m3_to_m8 x
m2_to_m7 x
A1_to_m7 x
m2_to_m8 x
A1_to_m8 
"""

data = '/Users/biggus/Documents/James/Data/Tu_miRNA/ResultsMiRNA_Targeting/2009_10_09run/2009_10_09runs.100bothCtrls.stats_1stdv.events.medFDRmeth.txt'
data = map(lambda l: l.strip('\n').split(' : '), open(data,'rU').readlines())

histOfField = 8  # <- field to make hist from in filtered data.

title = 'Hist of FDRs for combined (1stdv,class3)'
xAxisName = 'False Discovery Rate'
yAxisName = 'miRNAs'


def filterBy(data,index,srchTxt):
    """
    Returns list of data that contain <srchTxt> at position <index>.
    """
    rData = []
    for each in data:
        if each[index].find(srchTxt) != -1:
            rData.append(each)
    return rData


# write combination of filterBy iterations to get final data list

allTots = filterBy(data,1,'allPassedSeedsFor_')
lrgstFDR = max([eval(x[3]) for x in allTots])

filtData = filterBy(data,1,'allPassedSeedsFor_3')



consHist = [eval(x[3]) for x in filtData]
medHist  = [eval(x[4]) for x in filtData]


fig = plb.figure()
fig.suptitle(title)
ax = fig.add_subplot(111)

 

binSeq = [x*0.01 for x in range(int(round(100*lrgstFDR)+1))]
ax.hist(medHist,bins=binSeq,alpha=0.6,normed=0,cumulative=1,histtype='bar',label='medianFDR')
ax.hist(consHist,bins=binSeq,alpha=0.6,normed=0,cumulative=1,histtype='bar',label='medianFDR+sigma')
ax.set_xlabel(xAxisName)
ax.set_ylabel(yAxisName)
ax.text(0.8,0.8,
        'Total miRNAs = %s'\
        %(len(filtData)),
        bbox=dict(facecolor='gray',alpha=1),
        horizontalalignment='center',
        verticalalignment='center',
        transform= ax.transAxes)
ax.legend()
plb.show()


