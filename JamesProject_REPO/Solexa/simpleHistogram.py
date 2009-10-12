import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
#from pvalHist import normalize



inFile = '/Users/biggus/Downloads/s_7_sequence.seqCounts.txt'
inFile = open(inFile,'r')

counts = map(lambda l: int(l.strip()), inFile.readlines())

binNum = (max(counts)//50)+1

# the histogram of the data
n, bins, patches = plt.hist(counts, binNum, normed=1,facecolor='green', alpha=0.75, log=False)
#plt.semilogy(basey=10)
for i in range(0,9):
    print '%.4f' % ((bins[i+1]-bins[i]))

#bins = normalize(bins)

# add a 'best fit' line
#y = mlab.normpdf(bins, mu, sigma)
#l = plt.plot(bins)

plt.xlabel("Counts")
plt.ylabel('Sequences')
#plt.title('')
#plt.axis([320, 400, 0, 20])
plt.grid(True)

plt.show()

