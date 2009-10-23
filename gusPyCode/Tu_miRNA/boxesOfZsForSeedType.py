from matplotlib import pylab as plb
from gusPyCode.defs.miRNA_targeting import _seedModels
"""BoxPlots of 'z-scores' for each seedType"""



iFile = '/Users/biggus/Documents/James/Data/Tu_miRNA/ResultsMiRNA_Targeting/2009_08_31run.100ctrls.stats.txt'
iFile = map(lambda l: l.strip().split('\t'), open(iFile,'rU').readlines())

iFile.pop(0) # remove header

data = {}
# Initialize Data Dict
for seedType in _seedModels:
    data[seedType] = [[],[],[],[]]
    
    
# Extract Data
for line in iFile:
    data[line[1]][int(line[2])].append(float(line[10]))
    

fig = plb.figure()
ax1 = fig.add_subplot(4,1,1)
ax1.boxplot([data[x][0] for x in sorted(data.keys())], \
                    notch=0, sym='x', vert=1, whis=1.5,\
                    positions=None, widths=None,)
ax1.set_xticklabels([x for x in sorted(data.keys())])
ax1.set_ylabel('Absent in any.')
a
#ax1.grid(b=True)

ax2 = fig.add_subplot(4,1,2)
ax2.boxplot([data[x][1] for x in sorted(data.keys())], \
                    notch=0, sym='x', vert=1, whis=1.5,\
                    positions=None, widths=None,)
ax2.set_xticklabels([x for x in sorted(data.keys())])
ax2.set_ylabel('Present in 1')
#ax2.grid(b=True)

ax3 = fig.add_subplot(4,1,3)
ax3.boxplot([data[x][2] for x in sorted(data.keys())], \
                    notch=0, sym='x', vert=1, whis=1.5,\
                    positions=None, widths=None,)
ax3.set_xticklabels([x for x in sorted(data.keys())])
ax3.set_ylabel('Present in 2')
#ax3.grid(b=True)

ax4 = fig.add_subplot(4,1,4)
ax4.boxplot([data[x][3] for x in sorted(data.keys())], \
                    notch=0, sym='x', vert=1, whis=1.5,\
                    positions=None, widths=None,)
ax4.set_xticklabels([x for x in sorted(data.keys())])
ax4.set_ylabel('Present in 3')
#ax4.grid(b=True)

plb.show()
