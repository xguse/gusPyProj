from Bio.Blast import NCBIXML
from bioDefs import bestIdentOverLen

xmlOut = '/Users/biggus/Documents/James/Data/MicroArrayPrep/Aedes/Probes/ProbeBlastResults/Aedes.agilent.probes.nr.test.blastn_low.self.xml'


# load output
recs = NCBIXML.parse(open(xmlOut))

for rec in recs:
    print 'query_%s:' % (rec.query_id)
    for a in range(len(rec.alignments)):
        print '\tsubj_%s:' % (rec.alignments[a].hit_id)
        for h in range(len(rec.alignments[a].hsps)):
            print '\t%s\n\t%s\n\t%s\n' % \
                  (rec.alignments[a].hsps[h].query,rec.alignments[a].hsps[h].match,rec.alignments[a].hsps[h].sbjct)
        print '----------------------'
    print '========================================================================'
    None
    
    


None
