from TAMO import MotifMetrics, MotifTools
import TAMO
import math


totalSeqs   = MotifMetrics.ProbeSet('/Users/biggus/Documents/James/Data/2KB/2kb_Sequence/2kb_Anopheles/2KBupTSS_goodAffyAGAPsFastasOUT.masked.nr.fas')
clusterSeqs = TAMO.seq.Fasta.ids('/Users/biggus/Documents/James/Data/ClusterDefs/TC-Fastas/TC-2.fas')

motif       = MotifTools.Motif_from_text('WGATAAS')

e   = totalSeqs.Enrichment('WGATAAS',clusterSeqs)
ROC = 'ROC code not obtained yet'#totalSeqs.ROC_AUC('WGATAAS',clusterSeqs)


print "pVal =\t%s\nROC_AUC=\t%s" % (e,ROC)
