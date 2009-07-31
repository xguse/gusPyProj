from scipy import mat, transpose
import glob
import cPickle
from TAMO.MotifTools import Motif_from_counts,save_motifs

pklFile   = '/Users/biggus/Documents/James/Data/JasparMotifs/JASPAR_CORE_2008/jsprCore2008.tamo.pkl'
tmoFile   = '/Users/biggus/Documents/James/Data/JasparMotifs/JASPAR_CORE_2008/jsprCore2008.tamo.tmo'
jasparDir = '/Users/biggus/Documents/James/Data/JasparMotifs/JASPAR_CORE_2008/*.pfm'

nucBack = {'A':0.25,
           'C':0.25,
           'G':0.25,
           'T':0.25}

jasparFiles = glob.glob(jasparDir)


tamoMotifs = []

for jMat in jasparFiles:
    tempMat = map(lambda l: l.strip().split(), open(jMat, 'rU').readlines())
    
    ## eval() inteligently converts text numbers to int or float!
    #for i in range(len(tempMat)):
        #for j in range(len(tempMat[i])):
            #tempMat[i][j] = eval(tempMat[i][j]) 
    
    # transpose matrix
    tempMat = transpose(mat(tempMat)).tolist()
    
    for i in range(len(tempMat)):
        tempMat[i] = {'A':eval(tempMat[i][0]),
                      'C':eval(tempMat[i][1]),
                      'G':eval(tempMat[i][2]),
                      'T':eval(tempMat[i][3])}
    jasparTAMO_Motif = Motif_from_counts(tempMat,bg=nucBack)
    jasparTAMO_Motif.sourceFile = jMat.split('/')[-1]
    tamoMotifs.append(jasparTAMO_Motif)
    
# print to file
cPickle.dump(tamoMotifs,open(pklFile,'w'))
save_motifs(tamoMotifs,tmoFile,kmer_count=60)



print 'Done.'    

    
    