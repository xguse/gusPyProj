from gusPyCode.defs.bioDefs import cnvrtContigsInWig
from gusPyCode.defs.mosqData.AaCntgCvrt import supContigConvert

wigPath    = '/Users/biggus/Documents/James/Data/genomes/AaegL1/LBx_bowtie.minus_rRNA.wig'
outPath    = '/Users/biggus/Documents/James/Data/genomes/AaegL1/LBx_bowtie.minus_rRNA.cnvtd.wig'
cnvrtnDict = supContigConvert


cnvrtContigsInWig(wigPath,outPath,cnvrtnDict)