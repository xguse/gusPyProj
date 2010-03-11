import cPickle

# +++++ Definitions +++++
from gusPyCode.Tu_miRNA.mirTargetFuncGroups import getMasterTargetList 

def parseMicroArray(path2MicroArray):
    """"""
    maFile = open(path2MicroArray, 'rU')
    maList = []
    
    for line in maFile:
        maList.append(line.strip('\n').split('\t'))
        
    return maList
    
    
def sortGenesBy(maData,column):
    """"""
    maData.sort(key=lambda x: x[column])
    return maData
    
def plotKS(maData):
    """"""
    
    
# ++++++ hard coded variables +++++++
A1_m8_file = '/Users/biggus/Documents/James/Data/Tu_miRNA/ResultsMiRNA_Targeting/2009_11_11/2009_11_11.AGAP.stats.100.medFDRmeth.A1-m8passed.txt'
pklPath_Ca = '/Users/biggus/Documents/James/Data/Tu_miRNA/ResultsMiRNA_Targeting/2009_11_11/2009_11_11.AGAP.seedMatches.100mvCtrls.storeEvents.pkl'
# Set Event Class
useClass = 3

# -- Good miRNA A1-m8 list --
print "processing Good miRNA A1-m8 list..."
A1_m8s = sorted(map(lambda l: l.strip('\n'), open(A1_m8_file,'rU')))

# -- Process Events Pickle --
print "processing Events Pickle..."
data_Ca    = cPickle.load(open(pklPath_Ca,'rU'))
classConvert = {'II':2,'III':3}

ctrlNum = len(data_Ca[data_Ca.keys()[0]].ctrlEvents['A1_to_m8'])
masterGeneSets = [set(),initList(ctrlNum,set())]


getMasterTargetList(useClass,A1_m8s,masterGeneSets, data_Ca)
print 'deleting Events Pickle...'
del(data_Ca) # free-up some memory to work with!


