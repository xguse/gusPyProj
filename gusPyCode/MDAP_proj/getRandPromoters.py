import sys
import time
from TAMO.seq import Fasta

"""Takes a fastaFilePath and a fraction between 0 and 1.  Returns two fasta
files containing random sequences from fastaFilePath split randomly into files
of size 'fraction' and 1-'fraction'.  Files named as: fastaFilePath_fraction.Date_Time.fas"""

assert len(sys.argv[1:]) == 2, \
       'usage = %s fastaFilePath fraction<0 to 1>' % (sys.argv[0].split('/')[-1])
assert float(sys.argv[2]) <= 1 and float(sys.argv[2]) >= 0, \
       'usage = %s fastaFilePath fraction<0 to 1>' % (sys.argv[0].split('/')[-1])


filePath  = sys.argv[1]
frac      = float(sys.argv[2])
versionID = time.ctime().split(' ')
versionID = '%s%s_%s' % (versionID[1],versionID[2],versionID[3].replace(':','-'))


dict1,dict2 = Fasta.random_split(filePath,frac)

out1 = '%s_%s_%s.fas' % (filePath.split('/')[-1].rstrip('.fas'),frac,versionID)
out2 = '%s_%s_%s.fas' % (filePath.split('/')[-1].rstrip('.fas'),1-frac,versionID)

Fasta.write(dict1,out1,linelen=100)
Fasta.write(dict2,out2,linelen=100)
