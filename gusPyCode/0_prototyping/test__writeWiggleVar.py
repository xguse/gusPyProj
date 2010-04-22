from gusPyCode.defs.prototypes import _writeWiggleVar



#btPath  = '/home/dunnw/data/solexa/bowtie_out/LBx_bowtie.minus_rRNA.map_filtered.txt'
trkName = 'test_wigVarOut'
trkDesc = 'bowtie coverage after rRNA reads have been filtered out'
outName = '/Users/biggus/sandbox/testWriteWig/test_wigVarOut.wig'


aCnts = {'1':{3:4,
              4:3,
              5:5,
              6:7,
              10:10,
              11:9,}}



_writeWiggleVar(trkName, trkDesc, aCnts, outName)