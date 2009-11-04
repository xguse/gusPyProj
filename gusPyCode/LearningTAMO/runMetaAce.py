from TAMO.MD.AlignAce import MetaAce
from TAMO.seq import Fasta
from gusPyCode.defs.seqStats import calcStats
import glob
import sys
import pickle


assert 4 >= len(sys.argv[1:]) >= 3, 'usage: python %s fastaDir iterations runName [globPatern]' \
       % (sys.argv[0].split('/')[-1])

if len(sys.argv[1:]) > 3:
    globPattern = sys.argv[4]
else:
    globPattern = '*.fas'
    
fastaPaths = glob.glob(sys.argv[1]+globPattern)

gcBacks = []
for path in fastaPaths:
    stats = calcStats(path)
    gcBacks.append(stats['percentGC'])
    

metaAceObjs = []
for i in range(len(fastaPaths)):
    metaAceObjs.append(MetaAce(fastaPaths[i],iterations=int(sys.argv[2]),gcback=gcBacks[i]))
    
outFile = open('%s%s.pickle' % (sys.argv[1],sys.argv[3]),'w')
pickle.dump(tuple(metaAceObjs),outFile)
    



