from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Graphics import GenomeDiagram as gd
from reportlab.lib.units import cm



# ----- Fake up the SeqReq with SeqFeats -----
sFeats = [[[500,900,1],
           [34,40,1]],
          [[29,38,-1]]]


diag   = gd.Diagram('My Awesome Test Diagram')  # (1) open diag
for trx in sFeats:
    nTrk = diag.new_track(1,name='test',greytrack=1,greytrack_labels=1,scale=1) # (2) open track and set
    nSet = nTrk.new_set()
    for f in trx:
        nSet.add_feature(SeqFeature(FeatureLocation(f[0],f[1]),strand=f[2]),color=colors.blue, label=True) # (3) add features

## -- add trk1 features --
#fSet1.add_feature(sFeats[0],color=colors.blue, label=True)
#fSet1.add_feature(sFeats[1],color=colors.blue, label=True)
## -- add trk2 feats --
#fSet2.add_feature(sFeats[2],color=colors.red, label=True)

#diag.draw(format='linear', pagesize=(2*cm,19*cm), orientation='landscape', fragments=1, start=0, end=2000)
diag.draw(format='linear', fragments=1, start=0, end=2000)

diag.write("/Users/biggus/sandbox/TestBioGenomeDiagram.pdf", "pdf")

print 'Done.'
