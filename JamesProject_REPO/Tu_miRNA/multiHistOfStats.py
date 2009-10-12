from matplotlib import pylab as plb
"""plot multi hists for miRNA_targeting results to compare diff StDv metrics"""

iFile = ''
iFile = map(lambda l: l.strip().split('\t'), open(iFile,'rU').readlines())

iFile.pop(0) # remove header

groups = [[],[],[],[]]


for line in iFile:
    if   int(line[2]) == 0:
        groups[0].append(line)
    elif int(line[2]) == 1:
        groups[1].append(line)
    elif int(line[2]) == 2:
        groups[2].append(line)
    elif int(line[2]) == 3:
        groups[3].append(line)
        

fig = plb.figure()
ax1 = fig.subplot(1,4,1)
