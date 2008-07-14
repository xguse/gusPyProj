#!/sw/bin/python

from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
from Bio import SeqIO
import JamesDefs
import re
import string
from time import time


#--------- Script Specific Function Definitions ---------------------

def findAllMotifs(motifList, seqName, dictOfFastas, resultList):
    
    for motif in motifList:
        
        #  convert from IUPAC to regEx and search in forward direction
        fwd_RegExMotif = re.compile(JamesDefs.iupac2regex(motif), re.IGNORECASE)
        
        #  initiate result string for forward matches with name of AGAP and fwd IUPAC motif string
        fwd_MatchesStr = seqName+'\t'+motif+'\t'
        
        #  sequentially append each hit's coords to the end of fwd_MatchesStr
        for fwdMatcheObj in fwd_RegExMotif.finditer(dictOfFastas[seqName].seq.tostring()):
            
            ## must add 1 to the start pos due to computer numbers.............*
            fwd_MatchesStr = fwd_MatchesStr+'%s\t' % (str(fwdMatcheObj.start()+1)) 
            
        # add trailing newline for printing to file later
        fwd_MatchesStr = fwd_MatchesStr+'\n'
        
        # send fwd results to resultList
        resultList.append(fwd_MatchesStr)
        
        #  convert from IUPAC to regEx and search in reverse direction
        rev_RegExMotif = re.compile(JamesDefs.iupac2regex(JamesDefs.revComp(motif)), re.IGNORECASE)
        
        #  initiate result string for forward matches with name of AGAP and fwd IUPAC motif string
        rev_MatchesStr = seqName+'\t'+motif+'_rc\t'
        
        #  sequentially append each hit's coords to the end of rev_MatchesStr
        for revMatcheObj in rev_RegExMotif.finditer(dictOfFastas[seqName].seq.tostring()):
            
            ## must add 1 to the start pos due to computer numbers.............*
            rev_MatchesStr = rev_MatchesStr+'%s\t' % (str(revMatcheObj.start()+1)) 
            
        # add trailing newline for printing to file later
        rev_MatchesStr = rev_MatchesStr+'\n'
        
        # send fwd results to resultList
        resultList.append(rev_MatchesStr)

#--------------------------------------------------------------------



#========================= User Defined Variables =========================

#  InFiles:
motifList = map(string.strip, open('/Users/biggus/Documents/James/Data/2KB/2kb_Sequence/2kb_Combo/2Kb_AllMosquitoes/MosqMotifs/collapsed_upstream_ex-conserved_mosquito-motifs_nr.copRvCmp.txt', 'r'))

goodAGAPs = '/Users/biggus/Documents/James/Data/2KB/2kb_Sequence/2kb_Anopheles/2KBupTSS_goodAffyAGAPsFastasOUT.masked.nr.fas'

#  OutFile:
outFile   = '/Users/biggus/Documents/James/Data/2KB/2kb_Sequence/2kb_Combo/2Kb_AllMosquitoes/MosqMotifs/collapsed_upstream_ex-conserved_mosquito-motifs_nr.copRvCmp.map' 

#==========================================================================

t1 = time()
#  Create dict of fasta Objects from those in goodAGAPs
dictOfFastas = JamesDefs.fastaFileToBioSeqDict(goodAGAPs, Alphabet='IUPACAmbiguousDNA')


#  Create resultList to receive outPut
resultList = []

for seqName in dictOfFastas:
    
    findAllMotifs(motifList, seqName, dictOfFastas, resultList)
    lenResultList = len(resultList)

outFile= open(outFile, 'w')

outFile.writelines(resultList)

t2 = time()



print 'The operations took %.3f min to complete!' % ((float(t2)-t1)/60)
    