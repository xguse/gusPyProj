from TAMO.seq import Fasta
from gusPyCode.defs import bioDefs

miRNAFile = '/Users/biggus/Documents/James/Data/Tu_miRNA/miRNAs/miRBase/mature.aga.fa'
seedFile  = '/Users/biggus/Documents/James/Data/Tu_miRNA/miRNAs/miRBase/mature.aga.seeds.ctrl.fa'

oligoType = 'control' # 'match' or 'control'
assert oligoType == 'match' or 'control', 'oligoType MUST be only "match" or "control".'

# Load miRNA fastas into dict.
miRNAs = Fasta.file2dict(miRNAFile)

# Create new dict for seeds.
seeds = {}

# 1) Cycle through miRNA dict taking 7mers starting at pos 1 
#    and then pos2. Adapt key to reflect which. 
# 2) Convert to all uppers and convert U's to T's
# 3) If oligoType == 'match', rvcmp each 7mer and adapt key
#    to reflect which.
for miRNA in miRNAs:
    pos1_seed = miRNAs[miRNA][:7].upper().replace('U','T')
    pos2_seed = miRNAs[miRNA][1:8].upper().replace('U','T')


    if oligoType == 'match':
        seeds[miRNA+'_match_pos1'] = bioDefs.revComp(pos1_seed)
        seeds[miRNA+'_match_pos2'] = bioDefs.revComp(pos2_seed)
    else:
        seeds[miRNA+'_ctrl_pos1'] = pos1_seed
        seeds[miRNA+'_ctrl_pos2'] = pos2_seed
        
# Write out seed dict as fasta. 
Fasta.write(seeds,seedFile)

print "Done."