from JamesDefs import *


motif = 'WGATAR'

regExObj = iupac2fwdRevCmpTRegExObj(motif)


test = 'ATGCTAGCAGATAGTTACTACGTTATCTAAAAA'

if regExObj.search(test):
    print 'YES'
    
    
print 'done'

