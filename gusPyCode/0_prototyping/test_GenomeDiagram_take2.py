from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Graphics import GenomeDiagram as gd
from reportlab.lib.units import cm



# ----- Fake up the SeqReq with SeqFeats -----
sFeats = []
sFeats.append(SeqFeature(FeatureLocation(500,900),strand=1))
sFeats.append(SeqFeature(FeatureLocation(34,40),strand=1))
sFeats.append(SeqFeature(FeatureLocation(29,38),strand=-1))

diag   = gd.Diagram('My Awesome Test Diagram')
fTrak1 = diag.new_track(1,name='Track1',greytrack=1,greytrack_labels=1,greytrack_font_rotation=-30,scale=1)
fTrak2 = diag.new_track(1,name='Track2')

fSet1 = fTrak1.new_set()
fSet2 = fTrak2.new_set()

# -- add trk1 features --
fSet1.add_feature(sFeats[0],color=colors.blue, label=True)
fSet1.add_feature(sFeats[1],color=colors.blue, label=True)
# -- add trk2 feats --
fSet2.add_feature(sFeats[2],color=colors.red, label=True)

#diag.draw(format='linear', pagesize=(2*cm,19*cm), orientation='landscape', fragments=1, start=0, end=2000)
diag.draw(format='linear', fragments=1, start=0, end=2000)

diag.write("/Users/biggus/sandbox/TestBioGenomeDiagram.pdf", "pdf")

print 'Done.'
