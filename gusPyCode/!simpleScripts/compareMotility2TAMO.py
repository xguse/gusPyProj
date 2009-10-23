from TAMO.MotifTools import Motif
import motility

tM = Motif('WGATAR')
sites = tM.bogus_kmers()

tM = Motif(sites)
mM = motility.make_pwm(sites)


s = 'ATGCATGCTAGCGGCTGATAACGCTTATCATATGC'

mReults = mM.find(s,mM.max_score()*0.75,)