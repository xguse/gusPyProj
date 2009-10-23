from matplotlib import pylab as plb
from scipy import stats
import cPickle

inPkl = '/Users/biggus/Documents/James/Data/Tu_miRNA/ResultsMiRNA_Targeting/ctrlsComp_20090831_20090929.pkl'

# pickle in form:
#   [[data1_a],[data1_b],[data2_a],[data2_b] ...]
data = cPickle.load(open(inPkl,'rU'))

# Create new lists in data representing the combination of all subtypes
# THis is a custom HARD CODED area for now.
data.append(data[0]+data[2]+data[3])
data.append(data[1]+data[3]+data[5])

titOrder = iter(['Class I','Class II','Class III','All'])

for i in range(0,len(data),2):
    wxnStats = stats.wilcoxon(data[i],data[i+1])
    diffs    = [x[0]-x[1] for x in zip(data[i],data[i+1])]
    normStat = stats.normaltest(diffs)
    skew     = stats.skew(diffs)
    skewTest = stats.skewtest(diffs)
    
    fig = plb.figure()
    fig.suptitle('Paired Diffs for %s Controls' % (titOrder.next()))
    ax = fig.add_subplot(111)
    ax.hist(diffs,bins=30,normed=1)
    ax.set_xlabel('Difference Between Paired Controls (proSeed - Permutation)')
    ax.set_ylabel('Comparisons')
    ax.text(0.1,0.9,
            'Wilcoxon Signed-Rank Test:\nt-stat: %.2f\n2-tailed pVal: %.4G'\
            %(wxnStats[0],wxnStats[1]),
            bbox=dict(facecolor='gray',alpha = 0.5),
            horizontalalignment='left',
            verticalalignment='center',
            transform= ax.transAxes)
    ax.text(0.1,0.5,
            'd\'AugPear Norm Test:\nX^2:%.3f\n2-tailed pVal: %.4G\nSkew: %.4G\nSkew Test:\nz-Scr: %.4G\nz-Prob: %.4G'\
            %(normStat[0],normStat[1],skew,skewTest[0],skewTest[1]),
            bbox=dict(facecolor='gray',alpha = 0.5),
            horizontalalignment='left',
            verticalalignment='center',
            transform= ax.transAxes)
    
plb.show()

