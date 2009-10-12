from matplotlib import pylab as pb

orthoType = 1

# Figure Info
title='Class %s' % (orthoType)
xAx='Target sets predicted'
yAx='miRNAs'
save2=''


data = map(lambda l: l.strip().split('\t'),\
           open('/Users/biggus/Documents/James/Data/Tu_miRNA/ResultsMiRNA_Targeting/2009_08_31run.100ctrls.stats.events.all.cmltv.txt',\
                'rU'))




oData = []
for line in data:
    if (int(line[1]) == orthoType) and (int(line[2]) > 0):
        oData.append(len(line[4].split("),(")))
    #else:
        #oData.append(0)
        

fig = pb.figure()
fig.suptitle(title)
ax  = fig.add_subplot(111)
ax.set_xlabel(xAx)
ax.set_ylabel(yAx)
ax.hist(oData, bins=20, normed=0)

pb.show()