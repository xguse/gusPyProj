from TAMO import MotifTools 
from TAMO.seq import Fasta 
from TAMO.MotifMetrics import ProbeSet 
from TAMO.MD.AlignAce import AlignAce 
from TAMO.MD.MDscan import MDscan 
from TAMO.MD.Meme import Meme 
#from TAMO.DataSources import GO
from time import time


MDbg         = '/Users/biggus/Documents/James/Data/2KB/2kb_Sequence/2kb_Anopheles/2KBupTSS_goodAffyAGAPsFastasOUT.masked.nr.MD.bg'

outFile      = '/Users/biggus/Documents/James/Data/ClusterDefs/testTAMOmetrics.txt'


fastaPath    = '/Users/biggus/Documents/James/Data/ClusterDefs/TC-Fastas/TC-96.oneLine.fas'
clusterIDS   = Fasta.ids(fastaPath)
totalSeqs    = ProbeSet('/Users/biggus/Documents/James/Data/2KB/2kb_Sequence/2kb_Anopheles/2KBupTSS_goodAffyAGAPsFastasOUT.masked.nr.fas')  



print 'running MDscan...'
tMD_1 = time()
MDmotifs   = MDscan(fastaPath) #,bgfile=MDbg)
tMD_2 = time()
MD_time = tMD_2-tMD_1
print 'MDscan took %.5f sec == %.3f min.\nMDscan found %s motifs.' % (MD_time,MD_time/60.0, len(MDmotifs.motifs))




x=1