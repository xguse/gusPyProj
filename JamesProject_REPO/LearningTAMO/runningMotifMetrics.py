from TAMO import MotifTools 
from TAMO.seq import Fasta 
from TAMO import MotifMetrics
from TAMO.MD.AlignAce import AlignAce 
from TAMO.MD.MDscan import MDscan 
from TAMO.MD.Meme import Meme 
from TAMO import Clustering
#from TAMO.DataSources import GO
from time import time

TC8_path = '/Users/biggus/Documents/James/Data/ClusterDefs/TC-Fastas/TC-8.fas'
TC8_ids  = Fasta.ids(TC8_path)
TC8_seqs = Fasta.seqs(TC8_path)
allSeqs  = MotifMetrics.ProbeSet('/Users/biggus/Documents/James/Data/2KB/2kb_Sequence/2kb_Anopheles/2KBupTSS_goodAffyAGAPsFastasOUT.masked.nr.fas')

outFile  = '/Users/biggus/Documents/James/Data/ClusterDefs/TC-8_MotifMetrics.5-12.txt'

roughBestKmers = []

for i in range(6,10):
    imers = MotifMetrics.top_nmers_seqs(i,TC8_seqs)
    roughBestKmers.extend(imers)
    print '%s %smers found.' % (len(imers), i)
    
kmerMetrics = ['Kmer\thGeoPval\tBinomOverRep\n']
    
for kmer in roughBestKmers:
    hGeoPval = allSeqs.Enrichment(kmer, TC8_ids)
    binom   = allSeqs.overrep(kmer,TC8_ids)
    kmerMetrics.append('%s\t%s\t%s\n' % (kmer,hGeoPval,binom))
    
    
outFile = open(outFile,'w')
outFile.writelines(kmerMetrics)

print "Done."

