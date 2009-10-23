from MDAP_defs import transfacLike2tamoMotif
from TAMO.MotifMetrics import ProbeSet
import os


dirPath = '/Users/biggus/Documents/James/Data/ReClustering/PrelimData_Grant_Feb09/kMerPWMs/Clus2_kmerSearch.6-8mers.top10pcnt/'
probSet = ProbeSet('/Users/biggus/Documents/James/Data/2KB/2kb_Sequence/2kb_Anopheles/2KBupTSS_goodAffyAGAPsFastasOUT.masked.nr.fas')
coRegSeqs = '/Users/biggus/Documents/James/Data/ReClustering/kmedsPear33Clus50x_2/Clus2_247genes.genes.txt'
coRegSeqs = map(lambda l: l.strip(), open(coRegSeqs, 'rU').readlines())


os.chdir(dirPath)

files = os.listdir(dirPath)

pwms = {}
for f in files:
    if f.endswith(".pwm.txt"):
        motifName,tamoMotif = transfacLike2tamoMotif(f)
        pwms[motifName] = tamoMotif


# calc HyperGeo for all motifs
mNames = pwms.keys()
mNames.sort()

print '#threshold = 75% of max score\n#Name\tone_letter\tpVal\tfraction_of_gene_set'
for m in mNames:
    print '%s\t%s\t%s\t%s' % (m,
                              pwms[m].oneletter, 
                              probSet.p_value(pwms[m],coRegSeqs,factor=0.75), 
                              probSet.frac(pwms[m],coRegSeqs,factor=0.75))
print 'done'

