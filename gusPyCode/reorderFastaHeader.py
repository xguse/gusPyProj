from gusPyCode.defs import JamesDefs
import string


fastaFile =  open('/Users/biggus/Desktop/AnGambiaeProteinSeq.txt','rU').readlines()
outFile = open('/Users/biggus/Desktop/AnGambiaeProteinSeqByTranscript.fas','w')
outFile2 = open('/Users/biggus/Desktop/AnGambiaeProteinSeqByTranscript.delta.txt','w') 

delta = []

i = 0
for line in fastaFile:
    if line.startswith('>'):
        line = line.replace('>','')
        line = line.rstrip('\n')
        fields = line.split('|')
        fastaFile[i] = '>'+fields[4]+'|'+'|'.join(fields)+'\n'
        delta.append(fastaFile[i])
        i+=1
    else:
        i+=1
        
outFile.writelines(fastaFile)
outFile2.writelines(delta)

