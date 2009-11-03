from gusPyCode.MDAP_proj.MDAP_defs import transfacLike2tamoMotif
from TAMO.MotifMetrics import ProbeSet
from TAMO import MotifTools
import os


listPath = '/Users/biggus/Documents/James/Data/ReClustering/PrelimData_Grant_Feb09/kMerPWMs/Clus2_kmerSearch.6-8mers.top10pcnt/Group5.txt'
probSetReal = ProbeSet('/Users/biggus/Documents/James/Data/2KB/2kb_Sequence/2kb_Anopheles/2KBupTSS_goodAffyAGAPsFastasOUT.masked.nr.fas')
probSetShuf = ProbeSet('/Users/biggus/Documents/James/Data/2KB/2kb_Sequence/2kb_Anopheles/2KBupTSS_goodAffyAGAPsFastasOUT.masked.nr.seqsShuffled.fas')
coRegSeqs = '/Users/biggus/Documents/James/Data/ReClustering/kmedsPear33Clus50x_2/Clus2_247genes.genes.txt'
coRegSeqs = map(lambda l: l.strip(), open(coRegSeqs, 'rU').readlines())

# set real or shuffled genome
#probSet=probSetReal

motifs = map(lambda l: l.strip(), open(listPath, 'rU').readlines())

for m in range(len(motifs)):
    motifs[m] = MotifTools.Motif_from_text(motifs[m])


group = 'Group5'




# calc HyperGeo for all motifs


print '#threshold = 75% of max score\n#Group\tone_letter\tpVal\tfraction_of_gene_set\tpVal (shuffled)\tfraction_of_gene_set (shuffled)'
for m in motifs:
    print '%s\t%s\t%s\t%s\t%s\t%s' % (group,
                                      m.oneletter, 
                                      probSetReal.p_value(m,coRegSeqs,factor=0.75), 
                                      probSetReal.frac(m,coRegSeqs,factor=0.75),
                                      probSetShuf.p_value(m,coRegSeqs,factor=0.75), 
                                      probSetShuf.frac(m,coRegSeqs,factor=0.75))
print 'done'

