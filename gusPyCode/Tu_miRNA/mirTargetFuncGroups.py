"""mirTargetFuncGroups.py:
Compare the length of the union of gene target 
predictions made by all miRNA A1-m8 seed-types1 with the unions of the n control sets to produce 
a total number of gene targets belonging to each GO class (biological process). We will 
report those classes that produce results with FDR 0.25, the same as we used for the FDR_Ca 
"""


import csv
import cPickle
from pprint import pprint
import matplotlib
#matplotlib.use('TkAgg')
from matplotlib import pylab as pl


from gusPyCode.defs.JamesDefs import initList
from gusPyCode.defs import mathDefs


# --------- File Parsing ---------
useClass = 3

anoXcsv = '/Users/biggus/Documents/James/Data/GeneOntology/Ag-Sep-2009-Web.csv'
A1_m8_file = '/Users/biggus/Documents/James/Data/Tu_miRNA/ResultsMiRNA_Targeting/2009_11_11/2009_11_11.AGAP.stats.100.medFDRmeth.A1-m8passed.txt'
pklPath_Ca = '/Users/biggus/Documents/James/Data/Tu_miRNA/ResultsMiRNA_Targeting/2009_11_11/2009_11_11.AGAP.seedMatches.100mvCtrls.storeEvents.pkl'
figTit     = 'Biological Process Classification'
# -- AnoXcel --
print "\n\nprocessing AnoXcel..."
anoX_reader = csv.reader(open(anoXcsv, 'rU'))
anoXdata = {}
goInfo   = {}
for row in anoX_reader:
    if not row[0].startswith('AGAP'):continue # dont want header
    
    if row[0][:-3] in anoXdata:
        anoXdata[row[0][:-3]].append([row[0]]+row[137:142])
        goInfo[row[137:142][-2]] = row[137:142] # Create a library of goTermInfo
    else:
        anoXdata[row[0][:-3]] = [[row[0]]+row[137:142]]
        goInfo[row[137:142][-2]] = row[137:142] # Create a library of goTermInfo


# -- Good miRNA A1-m8 list --
print "processing Good miRNA A1-m8 list..."
A1_m8s = sorted(map(lambda l: l.strip('\n'), open(A1_m8_file,'rU')))

# -- Process Events Pickle --
print "processing Events Pickle..."
data_Ca    = cPickle.load(open(pklPath_Ca,'rU'))
classConvert = {'II':2,'III':3}

ctrlNum = len(data_Ca[data_Ca.keys()[0]].ctrlEvents['A1_to_m8'])
masterGeneSets = [set(),initList(ctrlNum,set())]

for miR in A1_m8s:
    name,cls = (miR.split(':')[0] , classConvert[miR.split(':')[1]])
    if cls == useClass:
        events = data_Ca[name].matchEvents['A1_to_m8'][cls]
        real_agaps=[]
        for geneSet in events:
            for gene in geneSet:
                if gene.startswith('AGAP'):
                    real_agaps.append(gene)
        masterGeneSets[0].update(real_agaps)
        for i in range(ctrlNum):
            events = data_Ca[name].ctrlEvents['A1_to_m8'][i][cls]
            ctrl_agaps=[]
            for geneSet in events:
                for gene in geneSet:
                    if gene.startswith('AGAP'):
                        ctrl_agaps.append(gene)
            masterGeneSets[1][i].update(ctrl_agaps)
print 'deleting Events Pickle...'
del(data_Ca) # free-up some memory to work with!



# --------- Defs ---------

def collectGeneNames(goTerm,masterGenes,anoXdict):
    """\tFor a GO term, query anoXdict for which genes in masterGenes(real & ctrls)
    are tagged with itself, and return a list([ set(real) , [sets(ctrl)] ]).
    
    RETURNS: a list([ set(realGenes) , [sets(ctrlGenes)] ])
    """
    rList = [set(),initList(ctrlNum,set())]
    
    notInAnoXcel = set()
    
    # -- Collect Reals --
    """Should Collect those gene Names not found in AnoXcel"""
    for gene in masterGenes[0]:
        try:
            anoXdata[gene]
        except KeyError:
            notInAnoXcel.add(gene)
            continue
        geneData = [x[-2] for x in anoXdata[gene]]
        if goTerm in geneData:
            rList[0].add(gene)
    
    # -- Collect Ctrls --
    """Should Collect those gene Names not found in AnoXcel"""
    for i in range(ctrlNum):
        for gene in masterGenes[1][i]:
            try:
                anoXdata[gene]
            except KeyError:
                notInAnoXcel.add(gene) 
                continue
            geneData = [x[-2] for x in anoXdata[gene]]
            if goTerm in geneData:
                rList[1][i].add(gene) 
                
    # Return
    return rList

def calcFDRStats(goDictEntry):
    realCount  = len(goDictEntry[0])
    ctrlCounts = [len(x) for x in goDictEntry[1]]
    
    FDRs = []
    for i in range(len(ctrlCounts)):
        ctrl = ctrlCounts[i]
        real = realCount
        
        if real == 0:
            return (None)
        else:
            fdr = float(ctrl)/real
            
        if fdr > 1:
            FDRs.append(1.0)
        else:
            FDRs.append(fdr)
            
    maths = mathDefs.stdDv(FDRs,'median')
    maths = list(maths)
    return [realCount]+maths

def labelFDRs(rects,FDRs):
    """Add FDR value to top of each Bar.
    """
    # attach FDR labels
    for i in range(len(rects)):
        height = rects[i].get_height()
        ax.text(rects[i].get_x()+rects[i].get_width()/2., 1.05*height, '%.3f'% (FDRs[i]),
                ha='center', va='bottom')
    


# --------- Main Body ---------
# -- Collect GO term Data --
print 'Collecting GO-term Data...'
goRsltDict = {}
for gTerm in goInfo:
    goRsltDict[gTerm] = collectGeneNames(gTerm,masterGeneSets,anoXdata)

#pprint(goRsltDict)
#exit('done')
# -- Calculate FDRs --
print 'Calculating FDRs...'
goFDRs = []
for gTerm in goRsltDict:
    fdrStats= calcFDRStats(goRsltDict[gTerm])
    if fdrStats:
        print fdrStats
        real,stDv,med = fdrStats
        if med+stDv <= 0.25: # enforce FDR filter
            goFDRs.append((gTerm,real,med+stDv))
goFDRs.sort(key=lambda x: x[1])


# -- Output Results --
print 'Outputing Results...'
fig = pl.figure()
ax = fig.add_subplot(111)

xVals   = range(len(goFDRs))
barVals = [x[1] for x in goFDRs]
rects = ax.bar(xVals,barVals)

ax.set_title(figTit)
ax.set_ylabel('Predicted targets')
ax.set_xticklabels( [x[0] for x in goFDRs] )

labelFDRs(rects,[x[2] for x in goFDRs])

pl.show()