"""
"""

import cPickle


from gusPyCode.defs.JamesDefs import initList,groupByField_silent,Bag
from gusPyCode.defs.statsDefs import cumHypergeoP

# ======== DEFINITIONS ==========
def collectMiRNA_totals(pathToFile):
    lines = open(pathToFile,'rU').readlines() 
    
    groupedLines = groupByField_silent(lines,0,sep=' : ')
    groupedLines.sort(key=lambda x: x[0][0])
    
    rDict = {}
    
    for miR in groupedLines:
        data = Bag({'name':miR[0][0],
                    'orthoTypes':[],
                    'AGAPs':initList(4,set())
                    })
        
        for line in miR:
            if line[1].startswith("allPassedSeedsFor_"):
                orthoType = int(line[1][-1])
                data.orthoTypes.append(orthoType)
                agaps = []
                for group in eval(line[-1]):
                    for gene in group:
                        if gene.startswith('AGAP'):
                            agaps.append(gene)
                data.AGAPs[orthoType].update(agaps)
        rDict[data.name]=data
    return rDict
                
def get_mCosmHits(miRNA,mCosmData):
    
    miRNA_data = []
    for line in mCosmData:
        if line[1] == miRNA:
            miRNA_data.append(line[11])
            
    miRNA_genes =  set()
    for Tx in miRNA_data:
        miRNA_genes.add(Tx[:-3])
        
    return miRNA_genes


def getAllSites(eventPkl,miRNA):
    total = []
    
    for seed in eventPkl[miRNA].matchEvents:
        total.extend(eventPkl[miRNA].matchEvents[seed][1])
        
    return set(total)
    

# ===== Load files =====
pklPath_Ca = '/Users/biggus/Documents/James/Data/Tu_miRNA/ResultsMiRNA_Targeting/2009_11_11/2009_11_11.AGAP.seedMatches.100mvCtrls.storeEvents.pkl'
data_Ca    = cPickle.load(open(pklPath_Ca,'rU'))

mCosm = map(lambda l: l.strip('\n').split('\t'), open('/Users/biggus/Documents/James/Data/Tu_miRNA/ResultsMiRNA_Targeting/Compare2MicroCosm/Anopheles/v5.txt.anopheles_gambiae.aga','rU'))

myResults = collectMiRNA_totals('/Users/biggus/Documents/James/Data/Tu_miRNA/ResultsMiRNA_Targeting/2009_11_11/2009_11_11.AGAP.stats.100.medFDRmeth.txt')




# ===== Main =====
classConvert = {2:'II',3:'III'}
outLines = []

for mir in sorted(myResults):
    for x in myResults[mir].orthoTypes:
        
        my_genes = myResults[mir].AGAPs[x]
        total = getAllSites(data_Ca,mir)
        
        n = get_mCosmHits(mir,mCosm)  ## positives in population (mCosm hits)
        i = n.intersection(my_genes)  ## positives in sample (intersection with mine)
        m = len(total)- len(n)        ## negatives in population (total - n)
        N = len(my_genes)             ## sample size (number in my set)

        cumPval = cumHypergeoP(len(n),len(i),m,N)
        
        outLines.append('%s\t%s\t%.4g\t%s\t%s\t%s\t%s' % \
                        (mir,
                         classConvert[x],
                         cumPval,
                         len(n),
                         len(i),
                         m,
                         N)
                        )

for l in outLines:
    print l
None



#for line in mCosm:
    #if line[1] == name:
        #mir1.append(line[11])
    



#mir1_genes = set()

#for Tx in mir1:
    #mir1_genes.add(Tx[:-3])


#mir1_genes
#len(mir1_genes)
#mir1_myGenes = set()

#for seed in data_Ca['aga-miR-1'].matchEvents:
    #mir1_myGenes.update(data_Ca['aga-miR-1'].matchEvents[seed][3])


#mir1_myGenes
#mir1_myGenes = list(mir1_myGenes)

#for i in range(len(mir1_myGenes)):
    #mir1_myGenes[i] =  mir1_myGenes[i][1]


#mir1_myGenes 
#mir1_myGenes = set(mir1_myGenes)
#mir1_genes.intersection(mir1_myGenes)
#mir1_genes.intersection(mir1_myGenes))
#mir1_genes.intersection(mir1_myGenes))
#mir1_genes.intersection(mir1_myGenes)
#len(mir1_genes.intersection(mir1_myGenes))/float(mir1_myGenes)
#len(mir1_genes.intersection(mir1_myGenes))/float(len(mir1_myGenes))
#len(mir1_myGenes)
#len(mir1_genes.intersection(mir1_myGenes))
#mir1_all = []
#32:
# for seed in data_Ca['aga-miR-1'].matchEvents:
#mir1_all.append(data_Ca['aga-miR-1'].matchEvents[seed][1])


#mir1_all
#len(mir1_all)
#mir1_all = []
#for seed in data_Ca['aga-miR-1'].matchEvents:
#37:
#mir1_all.extend(data_Ca['aga-miR-1'].matchEvents[seed][1])


#len(mir1_all)
#mir1_agaps = []
#40:
#each in mir1_all:
#if each.startswith('AGAP'):
    #mir1_agaps.append(each)
    

#len(mir1_agaps)
#mir1_agaps
#len(mir1_agaps)
#from scipy import  stats
#hypG = stats.hypergeom
##?hypG.cdf
##?hypG
#1 - hypG.cdf(15,2264,111)
#1 - hypG.cdf(15,2264,262,111)
#hypG.cdf(15,2264,262,111)
#from gusPyCode.defs import statsDefs
##?statsDefs.cumHypergeoP
#statsDefs.cumHypergeoP(262,15,2264-262,111)
#statsDefs.cumHypergeoP(111,15,2264-111,262)
#statsDefs.cumHypergeoP(262,15,2264-262,111)
#statsDefs.hypergeoP/
##?statsDefs.hypergeoP
#statsDefs.hypergeoP(262,15,2264-262,111)
#statsDefs.hypergeoP(111,15,2264-111,262)
#statsDefs.hypergeoP(111,15,2264-111,262L)
#_ip.magic("pwd ")
#_ip.magic("hist ")
#_ip.magic("hist 50")
#_ip.magic("hist -n 60")
#_ip.magic("hist 70")