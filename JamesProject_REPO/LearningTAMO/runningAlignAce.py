from TAMO import MotifTools 
from TAMO.seq import Fasta 
from TAMO.MotifMetrics import ProbeSet 
from TAMO.MD.AlignAce import AlignAce 
from TAMO.MD.MDscan import MDscan 
from TAMO.MD.Meme import Meme 
#from TAMO.DataSources import GO
from time import time

fastaPath    = '/Users/biggus/Documents/James/Data/ClusterDefs/TC-Fastas/TC-96.oneLine.fas'
clusterIDS   = Fasta.ids(fastaPath)
totalSeqs    = ProbeSet(fastaPath)  # !! this is wrong should proly be goodAffys

MDbg         = '/Users/biggus/Documents/James/Data/2KB/2kb_Sequence/2kb_Anopheles/2KBupTSS_goodAffyAGAPsFastasOUT.masked.nr.MD.bg'

outFile      = '/Users/biggus/Documents/James/Data/ClusterDefs/testTAMOmetrics.txt'

#theAce = AlignAce(fastaPath,width=10)

print 'running MDscan...'
tMD_1 = time()
MDmotifs   = MDscan(fastaPath) #,bgfile=MDbg)
tMD_2 = time()
MD_time = tMD_2-tMD_1
print 'MDscan took %.5f sec == %.3f min.\nMDscan found %s motifs.' % (MD_time,MD_time/60.0, len(MDmotifs.motifs))

print 'running MEME...'
tMeme_1 = time()
memeMotifs = Meme(fastaPath)
tMeme_2 = time()
Meme_time = tMeme_2-tMeme_1
print 'Meme took %.5f sec == %.3f min.\nMeme found %s motifs.' % (Meme_time, Meme_time/60.0, len(memeMotifs.motifs))


print 'running metrics on MDscan motifs...'
tMetMD_1 = time()
print 'running "church" metric...'
for m in MDmotifs.motifs:
    m.church = totalSeqs.church(m, clusterIDS)
print 'running HyperG metric...'
for m in MDmotifs.motifs:
    m.hyperG = totalSeqs.Enrichment(m, clusterIDS)
tMetMD_2 = time()
MDmetTime = tMetMD_2-tMetMD_1
print 'Metrics on MDscan motifs took %.5f sec == %.3f min.' % (MDmetTime,MDmetTime/60.0)


    
print 'running metrics on MEME motifs...'
tMetMeme_1 = time()
print 'running "church" metric...'
for i in range(0, len(memeMotifs.motifs)):
    memeMotifs.motifs[i].church = totalSeqs.church(memeMotifs.motifs[i], clusterIDS)
    x=1
print 'running HyperG metric...'
for m in memeMotifs.motifs:
    m.hyperG = totalSeqs.Enrichment(m, clusterIDS)
tMetMeme_2 = time()
MemeMetTime = tMetMeme_2-tMetMeme_1
print 'Metrics on Meme motifs took %.5f sec == %.3f min.' % (MemeMetTime,MemeMetTime/60.0)


outFile = open(outFile,'w')
outFile.write('# MDscan output:\nmotif\tchurch\thyperGeo\n')
for m in MDmotifs.motifs:
    outFile.write('%s\t%s\t%s\n' % (m.motif,m.church,m.hyperG))
    
outFile.write('# Meme output:\nmotif\tchurch\thyperGeo\n')
for m in memeMotifs.motifs:
    outFile.write('%s\t%s\t%s\n' % (m,m.church,m.hyperG))

x=1