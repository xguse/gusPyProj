import re


#========================= User Defined Variables =========================
fastaFile   = ''
peptideFile = ''
outFile     = ''

#==========================================================================

#--------- Script Specific Function Definitions ---------------------
def upperIt(matchObj):
    return matchObj.group(0).upper()

#--------------------------------------------------------------------

# (0) open and initiate data lists
fastaFile = map(lambda line: line.strip(), open(fastaFile, 'rU').readlines())
peptideFile = map(lambda line: line.strip(), open(peptideFile, 'rU').readlines())




# (1) convert file to lowercase
for i in range(0,len(fastaFile)):
    
    if fastaFile[i][0] == '>':
        continue
    else:
        fastaFile[i] = fastaFile[i].lower()
del i


# (2) modify peptide data list to accomadate sequence counts by converting to list of lists [[pep1,count],[pep2,count]]
for i in range(0, len(peptideFile)):
    peptideFile[i] = [peptideFile[i], 0]
    
    

# (3)loop through fasta file converting matches to UPPERcase
for i in range(0,len(inFile)):
    
    if fastaFile[i][0] == '>':
        continue
    else:
        for peptide in peptideFile:
            pepEx = re.compile(peptide[0], re.IGNORECASE)
            
            fastaFile[i] = re.sub(pepEx, upperIt, fastaFile[i])