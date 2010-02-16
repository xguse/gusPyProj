from weblogolib import *
#import pickle

fin = open('/Users/biggus/SourceCodeEtc/weblogo-3.0/tests/data/cap.fas')
seqs = read_seq_data(fin)
data = LogoData.from_seqs(seqs)
#pickle.dump(data,open('/Users/biggus/SourceCodeEtc/weblogo-3.0/tests/out/fromFasta.data.pkl','w'))
#exit('Just Dumped Fasta Pkl')
options = LogoOptions()
options.title = 'A Logo Title'
format = LogoFormat(data, options)
fout = open('/Users/biggus/SourceCodeEtc/weblogo-3.0/tests/out/cap.xmpl.pdf', 'w')
pdf_formatter( data, format, fout)
None
print 'done'