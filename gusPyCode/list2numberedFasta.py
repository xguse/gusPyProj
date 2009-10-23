#========================= User Defined Variables =========================
mofifInFile        = '/Users/biggus/Documents/James/Data/motifIdentitySearch/CulexAnopheles/CulexAnophUp8mer_gte3.doug.motifs'

motifFastaFile     = '/Users/biggus/Documents/James/Data/motifIdentitySearch/CulexAnopheles/CulexAnophUp8mer_gte3.doug.motifs.fas'

#==========================================================================

mofifInFile = map(lambda line: line.strip(), open(mofifInFile, 'rU').readlines())

fastaList = []

commonHeaderInfo = 'CulexAnoph8doug'

c = 1
for each in mofifInFile:
    
    record = '>%i:%s\n%s\n' % (c, commonHeaderInfo, each)
    fastaList.append(record)
    
    c+=1
    
    
motifFastaFile = open(motifFastaFile, 'w')
motifFastaFile.writelines(fastaList)

print '%i motifs written as fasta records.' % len(fastaList)