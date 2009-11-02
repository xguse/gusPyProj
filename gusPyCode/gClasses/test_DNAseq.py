from gusPyCode.gClasses.gSeqClasses import DNAseq
from Bio.Seq import Seq

q = Seq('kjhas')

d = DNAseq('''atgatgaagatcacgatcgccgccctcgctctgctggccgtcgccgttgctgcatacgaa
catgaggactatcattcgcatcccagctacaagttcgagtatggagtaaaggatcctcac
accggagaccacaagagccagtgggaacatcgggacggagatgtcgtcaagggagcatat
acccttcacgaggctgacggaactgagcgagtagttgagtactcgtccgacaagcacaac
ggcttccaggctcacgtcaagcgagtgggtcatgctcaccacccggaagtgtatggacac
catgaggctggacactcgtacggccatggacacggtcatggtcatggcagcagctacgcc
aacggcaacctgtaccagcaccatcactaa
''','testSeq1')

print 'repr:\n%s' % (d)
print 'str(d):\n%s' % (str(d))
print 'len(d):\n%s' % (len(d))
print 'toString:\n%s' % (d.toString())
print 'toFasta:\n%s'  % (d.toFasta())

# ResFreqs and stuff
print 'Freq of "A":\n%s' % (d.getResFreq('A'))
print 'FreqTable:\n'
frqTab = d.getResFreqTable()
for res in frqTab:
    print res


#slice
sd = d.slice(3,9)
print 'slice:\n%s' % (sd.toString())

# rvCmp
rd = d.revCmp()
print 'revCmp:\n%s' % (rd.toString())

# __add__
add_d = d+sd
print 'add:\n%s + %s = %s' % (d.toString(),sd.toString(),add_d.toString())

# codons
print 'codons1:\n%s'  % (d.codons())
print 'codons2:\n%s'  % (d.codons(2))
print 'codons3:\n%s'  % (d.codons(3))
print 'codons_1:\n%s' % (d.codons(-1))
print 'codons_2:\n%s' % (d.codons(-2))
print 'codons_3:\n%s' % (d.codons(-3))

# translate 
print 'translate testSeq1 in 1:\n%s'  % (d.translate().toString())
print 'translate testSeq1 in 2:\n%s'  % (d.translate(2).toString())
print 'translate testSeq1 in 3:\n%s'  % (d.translate(3).toString())
print 'translate testSeq1 in -1:\n%s' % (d.translate(-1).toString())
print 'translate testSeq1 in -2:\n%s' % (d.translate(-2).toString())
print 'translate testSeq1 in -3:\n%s' % (d.translate(-3).toString())
