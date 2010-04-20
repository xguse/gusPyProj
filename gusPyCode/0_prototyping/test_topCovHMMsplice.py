from gusPyCode.defs.bioDefs import topCovHMMsplices
from gusPyCode.defs.HMMsplice_utils import readJunctionsFromBed

jncPath = '/Users/biggus/Documents/James/Data/Solexa/aedes/hmmSplicer/finalResults/LB/junction.renamed.bed'

jDict = readJunctionsFromBed(jncPath,saveWholeLine=1)

topDict = topCovHMMsplices(jDict,frac=0.01)

