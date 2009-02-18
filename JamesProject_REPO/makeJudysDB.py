from TAMO.seq import Fasta

fastaFile = '/Users/biggus/Documents/James/AedesPeptides/Aedes_aegypti.AaegL1.50.pep.all.fa'
outFile   = '/Users/biggus/Documents/James/AedesPeptides/Aedes_aegypti.AaegL1.50.pep.all.reformatted.fa'

fDict = Fasta.load(fastaFile, lambda x: x)

newDict = {}

for n,s in fDict.items():
    n = [n[:13]+'_Ens',n[13:]]
    n = '|'.join(n)
    n = '|'+n
    
    newDict[n] = s
    
outFile = open(outFile, 'w')
newDict_keys = newDict.keys()
newDict_keys.sort()

for key in newDict_keys:
    entry = '>%s\n%s\n' % (key,newDict[key])
    outFile.write(entry)
    
print 'Done.'    