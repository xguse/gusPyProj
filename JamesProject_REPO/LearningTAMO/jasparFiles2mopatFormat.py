from scipy import mat, transpose
import glob

outFile   = '/Users/biggus/Documents/James/Data/JasparMotifs/JASPAR_CORE_2008/jsprCore2008.mopat.motifs'
jasparDir = '/Users/biggus/Documents/James/Data/JasparMotifs/JASPAR_CORE_2008/*.pfm'


jasparFiles = glob.glob(jasparDir)


jasparPfmDict = {}

for jMat in jasparFiles:
    tempMat = map(lambda l: l.strip().split(), open(jMat, 'rU').readlines())
    
    ## eval() inteligently converts text numbers to int or float!
    #for i in range(len(tempMat)):
        #for j in range(len(tempMat[i])):
            #tempMat[i][j] = eval(tempMat[i][j]) 
    
    # transpose matrix
    tempMat = transpose(mat(tempMat)).tolist()
    
    jasparPfmDict[jMat.split('/')[-1].replace('.pfm','')] = tempMat
    
# print to file
keys = sorted(jasparPfmDict.keys())

outFile = open(outFile, 'w')
for k in keys:
    toOut = '>%s\t\n' % (k)
    for i in range(len(jasparPfmDict[k])):
        toOut += '%s\n' % (' '.join(jasparPfmDict[k][i]))
    outFile.write(toOut)
    
outFile.close()

print 'Done.'    

    
    