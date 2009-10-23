from TAMO import MotifMetrics, MotifTools
from TAMO.seq import FakeFasta

fastaPath   = '/Users/biggus/Documents/James/Data/2KB/2kb_Sequence/2kb_Anopheles/2KBupTSS_goodAffyAGAPsFastasOUT.masked.nr.fas'
totalSeqs   = MotifMetrics.ProbeSet(fastaPath)

r1 = FakeFasta.random_seqs(3,fastaPath)
r2 = FakeFasta.random_seqs(3,fastaPath)

x=1
