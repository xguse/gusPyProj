from defs_moduleByTableLookUp import makeFwdAndRevCompRegExObj
from time import time

#--------- Script Specific Function Definitions ---------------------
def findMatchesInFastaDB(motif, DB):
    t1 = time()
    
    motifRegEx = makeFwdAndRevCompRegExObj(motif)
    
    hitList = []
    
    for i in range(0, len(DB)):
        if DB[i][0] == '>':
            continue
        else:
            if motifRegEx.search(DB[i]):
                hitList.append(DB[i-1])
                
    t2 = time()
    print '%s was found in %d records and the search took %.4f seconds' % (motif, len(hitList), t2-t1)  # %.4f
    return hitList
    

#--------------------------------------------------------------------



#========================= User Defined Variables =========================
mofifInFile = '/Users/biggus/Documents/James/Data/DougEmail/shuffled_Upstream_mosquito-conserved_8-mers.txt'

fastaDB     = '/Users/biggus/Documents/James/Data/RedFly/redFlyDoug.fas'

#outFile     = '/Users/biggus/Documents/James/Data/RedFly/shuffled_Upstream_mosquito-conserved_8-mers.redFlyCount.txt'
#==========================================================================

mofifInFile = map(lambda line: line.strip(), open(mofifInFile, 'rU').readlines())

fastaDB = map(lambda line: line.strip(), open(fastaDB, 'rU').readlines())




for each in mofifInFile:
    listOfMatchingFastaHeads = findMatchesInFastaDB(each, fastaDB)
    for each in listOfMatchingFastaHeads:
        print each
    print '\n'
    
    

