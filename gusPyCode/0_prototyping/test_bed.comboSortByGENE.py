from gusPyCode.defs.bed import comboSortByGENE
from gusPyCode.defs.bed import getMultiGene



tblA = '/Users/biggus/Documents/James/Data/Solexa/aedes/hmmSplicer/finalResults/LB/LBjunction.edges.NOT.LSedges.AND.aaegypti.Tx-Ensembl.bed' 
tblB = '/Users/biggus/Documents/James/Data/Solexa/aedes/hmmSplicer/finalResults/LS/LSjunction.edges.NOT.LBedges.AND.aaegypti.Tx-Ensembl.bed'
oPath = '/Users/biggus/Documents/James/Data/Solexa/aedes/hmmSplicer/finalResults/testSrtByGene.txt'
colA = 8
colB = 8

combo = comboSortByGENE(tblA,tblB,colA,colB,oPath)

multiTx = getMultiGene(combo,1)
mOut = open('/Users/biggus/Documents/James/Data/Solexa/aedes/hmmSplicer/finalResults/testMultiGroup.txt' , 'w')
for l in multiTx:
	mOut.write('%s\n' % ('\t'.join(l)))
mOut.close()