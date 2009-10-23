from TAMO.seq import Fasta

fastaFile = '/Users/biggus/Documents/James/AedesPeptides/aaegypti.PEPTIDES-AaegL1.1.reformatted.fa'

fDict = Fasta.load(fastaFile, lambda x: x.split('|')[1])

nrNames = []
rNames  = []
for each in fDict.keys():
    if each not in nrNames:
        nrNames.append(each)
    else:
        if each not in rNames:
            rNames.append(each)
            
print '%s names were repeated at least once.' % (len(rNames))

x=1
