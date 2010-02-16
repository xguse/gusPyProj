import weblogolib as wll
import pickle

seqs = ['attcgtgatagctgtcgtaaag','ttttgttacctgcctctaactt','aagtgtgacgccgtgcaaataa','tgccgtgattatagacactttt','atttgcgatgcgtcgcgcattt','taatgagattcagatcacatat','taatgtgacgtcctttgcatac','gaaggcgacctgggtcatgctg','aggtgttaaattgatcacgttt','cgatgcgaggcggatcgaaaaa','aaattcaatattcatcacactt','agatgtgagccagctcaccata','agatgtgattagattattattc','aattgtgatgtgtatcgaagtg','ttatttgaaccagatcgcatta','aaatgtaagctgtgccacgttt','aagtgtgacatggaataaatta','ttgtttgatttcgcgcatattc','aaacgtgatttcatgcgtcatt','atgtgtgcggcaattcacattt','taatgttatacatatcactcta','ttttatgacgaggcacacacat','aagttcgatatttctcgttttt','ttttgcgatcaaaataacactt','aaacgtgatcaacccctcaatt','taatgtgagttagctcactcat','aattgtgagcggataacaattt','ttgtgtgatctctgttacagaa','TAAtgtggagatgcgcacaTAA','TTTtgcaagcaacatcacgAAA','GACctcggtttagttcacaGAA','aattgtgacacagtgcaaattc','aaccgtgctcccactcgcagtc','TCTTGTGATTCAGATCACAAAG','ttttgtgagttttgtcaccaaa','aaatgttatccacatcacaatt','ttatttgccacaggtaacaaaa','atgcctgacggagttcacactt','taacgtgatcatatcaacagaa','Ttttgtggcctgcttcaaactt','ttttatgatttggttcaattct','aattgtgaacatcatcacgttc','ttttgtgatctgtttaaatgtt','agaggtgattttgatcacggaa','atttgtgagtggtcgcacatat','gattgtgattcgattcacattt','gtgtgtaaacgtgaacgcaatc','aactgtgaaacgaaacatattt','TCTTGTGATGTGGTTAACCAAT']

#fIn = open('/Users/biggus/SourceCodeEtc/weblogo-3.0/tests/data/cap.fas')
#seqs = wll.read_seq_data(fIn)


for i in range(len(seqs)):
    seqs[i] = wll.Seq(seqs[i])
    
seqs = wll.SeqList(seqs)
seqs.alphabet = wll.which_alphabet(seqs)

data = wll.LogoData.from_seqs(seqs)
#pickle.dump(data,open('/Users/biggus/SourceCodeEtc/weblogo-3.0/tests/out/fromList.data.pkl','w'))
#exit('Just Dumped List Pkl')
options = wll.LogoOptions()
options.title = 'A Logo Title'
format = wll.LogoFormat(data, options)
fout = open('/Users/biggus/SourceCodeEtc/weblogo-3.0/tests/out/cap2.xmpl.pdf', 'w')
wll.pdf_formatter( data, format, fout)
None
print 'done'