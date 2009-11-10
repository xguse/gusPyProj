from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
from Bio import SeqIO
from gusPyCode.defs import JamesDefs
import re
import string
from gusPyCode.defs.defs_moduleByTableLookUp import findAllMotifsInAGAP
from time import time

boundarySeqs          = '/Users/biggus/Documents/MBGB/Rotations/James/Data/testData4Mapper/MGandUpAt24.fas'
#  /Users/biggus/Documents/MBGB/Rotations/James/Data/2KB/2kb_Sequence/2kb_Anopheles/2kb_anophelesProcessed/2KB_anophelesUpstream/2KBupTSS_goodAffyAGAPsFastasOUT.masked.nr.fas  
motifList             = map(string.strip, open('/Users/biggus/Documents/MBGB/Rotations/James/Data/DegenMotifs/aedesAnopheles7mer.TSS.nr.motif', 'r'))
outFile               = open('/Users/biggus/Documents/MBGB/Rotations/James/Data/MotifMaps/matrixOfMGandUpAt24.txt','w')

#  Populate a Dict with Seq objs for Anopheles boundary seqs
#  What follows directly is a klugde to get my seqDict vals to have the IUPAC ambiguous alphabet
boundarySeqs = list(SeqIO.parse(open(boundarySeqs, "rU"), "fasta"))
for record in boundarySeqs :
    record.seq.alphabet = IUPACAmbiguousDNA

boundarySeqs = SeqIO.to_dict(boundarySeqs, key_function = lambda rec : rec.description.split()[0])

matrixByAGAP= {}

t1 = time()
for AGAP in boundarySeqs:
    tIN1= time()
    findAllMotifsInAGAP(motifList, AGAP, boundarySeqs, matrixByAGAP)
    x=0
    tIN2 = time()
    print 'findAllMotifsInAGAP for %s took %s seconds.' % (AGAP, str(tIN2-tIN1))

    
t2 = time()

print 'the matrix took %s seconds to build.' % (str(t2-t1))


outList = []
outList.append('\t'+'\t'.join(motifList)+'\n')

for AGAP in matrixByAGAP:
    lineString = AGAP
    for x in matrixByAGAP[AGAP]:
        lineString+='\t'+str(x)
    outList.append(lineString+'\n')
    x=0
    
outFile.writelines(outList)