from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Graphics import GenomeDiagram
from reportlab.lib.units import cm

gdd = GenomeDiagram.Diagram('Test Diagram')
gdt_features = gdd.new_track(1, greytrack=False)
gds_features = gdt_features.new_set()

#Add three features to show the strand options,
feature = SeqFeature(FeatureLocation(25, 125), strand=+1)
gds_features.add_feature(feature, name="Forward", label=True)
feature = SeqFeature(FeatureLocation(150, 250), strand=None)
gds_features.add_feature(feature, name="Standless", label=True)
feature = SeqFeature(FeatureLocation(275, 375), strand=-1)
gds_features.add_feature(feature, name="Reverse", label=True)
#-------------
#Large font, parallel with the track
gds_features.add_feature(feature, label=True, color="green",
                           label_size=25, label_angle=0)

#Very small font, perpendicular to the track (towards it)
gds_features.add_feature(feature, label=True, color="purple",
                           label_position="end",
                           label_size=4, label_angle=90)

#Small font, perpendicular to the track (away from it)
gds_features.add_feature(feature, label=True, color="blue",
                           label_position="middle",
                           label_size=6, label_angle=-90)
#-------------


gdd.draw(format='linear', pagesize=(15*cm,4*cm), fragments=1,
         start=0, end=400)
gdd.write("/Users/biggus/sandbox/GD_labels_default.pdf", "pdf")
