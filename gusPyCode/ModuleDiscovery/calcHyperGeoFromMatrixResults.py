import sys
from gusPyCode.defs.doug_hypergeometric import hyperGeoPvalue
from time import time

#========================= User Defined Variables =========================
paramsFile  = sys.argv[1]

outFile     = open(sys.argv[2],'w')
#===============================================

hyperGeoParams = map(lambda line : line.strip(), open(paramsFile, 'rU').readlines())

t1 = time()
c = 0
while c < len(hyperGeoParams):
    p1 = time()
    #  explode string to list of params
    params = hyperGeoParams[c].split('\t')
    
    #  calc hyperGeo pVal
    pVal = hyperGeoPvalue(int(params[1]),int(params[2]),int(params[3]),int(params[4]))
    
    outFile.write(hyperGeoParams[c]+'\t'+str(pVal)+'\n')
    outFile.flush()
    c+=1  
    p2 = time()
    print 'pval %s: %.6f  sec: %.8f' % (c,pVal,p2-p1)
t2 = time()

print 'Calculateing the hyperGeo pVals took %.2f min.' % ((t2-t1)/60.0)


