from Bio import SeqIO


filename = "/Users/biggus/Documents/MBGB/Rotations/James/Data/Sequence/Anopheles/testFasta1.fas"
d = SeqIO.to_dict(SeqIO.parse(open(filename, "rU"), 'fasta'),
    key_function = lambda rec : rec.description.split()[0])
print len(d)
print d.keys()[0:10]
key = d.keys()[0]
print d[key]
print d['AGAP000028'].seq[0:]