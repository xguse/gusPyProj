from gusPyCode.defs.HMMsplice_utils import *
import pyroc as roc

def toROCdata(scoredDict):
	rocData = []
	for k in scoredDict:
		rocData.append((scoredDict[k][1],scoredDict[k][0],k))
	return rocData

anno = '/Users/biggus/Documents/James/Data/genomes/AaegL1/aaegypti.Tx-Ensembl.bed'
mult = '/Users/biggus/Documents/James/Data/Solexa/aedes/hmmSplicer/finalResults/Lx_unfiltered/LX.gtag.collapsed.multi.bed'
sngl = '/Users/biggus/Documents/James/Data/Solexa/aedes/hmmSplicer/finalResults/Lx_unfiltered/LX.gtag.collapsed.sngl.bed'
mult = toROCdata(generateROCcurve(mult, anno, stepSize=10, wiggle=3)[1])
sngl = toROCdata(generateROCcurve(sngl, anno, stepSize=10, wiggle=3)[1])
mROC = roc.ROCData(mult,linestyle='r-')
sROC = roc.ROCData(sngl,linestyle='b-')
roc.plot_multiple_roc([mROC,sROC],title='Multiples vs Singles',labels=['Multiples','Singles'], include_baseline=1, equal_aspect=True)

mROC.auc()
sROC.auc()
