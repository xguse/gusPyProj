from TAMO import MotifTools

motifs = ['TGATACA',
          'TGATAAA',
          'TGACAAA',
          'AGATACA',
          'AGATACG',]
for i in range(len(motifs)):
    motifs[i] = MotifTools.Motif_from_text(motifs[i])

weights = [3,
           7,
           2,
           2,
           1,]

noWeight = MotifTools.sum(motifs)
withWeights = MotifTools.sum(motifs, weights)

print "nW:\n%s\n" % (noWeight.printlogo())
print "wW:\n%s" % (withWeights.printlogo())
x=1
