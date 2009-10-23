from gusPyCode.MDOSX_proj.MDOSX_defs.MDOSX_classes import OrthoGroup
from TAMO import MotifTools


orthoDict = 1
orthoDict = {'AGAP1':'AAAAAGATATTTTT',
             'AAEL1':'AAAAAGATATTTTT',
             'CPJI1':'AAAAAGATTTTTTT'}

oGroup = OrthoGroup(orthoDict)

motif = MotifTools.Motif_from_text('GATA')

oGroup.searchMotif(motif,scoreFactor='0.1')

None


