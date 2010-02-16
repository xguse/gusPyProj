from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Graphics import GenomeDiagram
from reportlab.lib.units import cm

gdd = GenomeDiagram.Diagram('Test Diagram')
tk1_features = gdd.new_track(1, greytrack=False)
tk1_st1_features = tk1_features.new_set()

tk2_features = gdd.new_track(2, greytrack=False)
tk2_st1_features = tk2_features.new_set()

tk1_st1_features = tk1_features.new_set()
tk1_st1_features.add_feature(SeqFeature(FeatureLocation(104, 117),color='green', strand=None))
tk1_st1_features.add_feature(SeqFeature(FeatureLocation(25, 125), color='red',strand=None))
tk1_st1_features.add_feature(SeqFeature(FeatureLocation(1666, 1673),color='blue', strand=None))



tk2_st1_features = tk2_features.new_set()
tk2_st1_features.add_feature(SeqFeature(FeatureLocation(150, 157), color='green',strand=None))
tk2_st1_features.add_feature(SeqFeature(FeatureLocation(10, 17), color='red',strand=None))
tk2_st1_features.add_feature(SeqFeature(FeatureLocation(1000, 1010), color='blue',strand=None))




gdd.draw(format='linear', pagesize=(2*cm,19*cm), orientation='landscape', fragments=1, start=0, end=2000)

gdd.write("/Users/biggus/sandbox/TestBioGenomeDiagram.pdf", "pdf")

print 'Done.'
