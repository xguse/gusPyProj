"""mirTargetFuncGroups.py:
Compare the length of the union of gene target 
predictions made by all miRNA A1-m8 seed-types1 with the unions of the n control sets to produce 
a total number of gene targets belonging to each GO class (biological process). We will 
report those classes that produce results with FDR 0.25, the same as we used for the FDR_Ca 
"""


import csv
import cPickle
from pprint import pprint
import matplotlib as mpl
mpl.use('WXAgg')
from matplotlib import pylab as pl


from gusPyCode.defs.JamesDefs import initList,Bag
from gusPyCode.defs import mathDefs


# --------- Defs ---------

def getMasterTargetList(oTyp2use,A1_m8List,masterDest, mirResultsPkl):
    """populate masterDest with the cumulative set of genes predicted to be
    targeted by the A1-m8 seed-types listed in A1_m8List.
    """
    # for seed in seedList,
    #    Add real_targets to cumulative real_list
    #    Add ctrl_targets_i to ctrl_list_i
    classConvert = {'II':2,'III':3}
    
    for mir in A1_m8List:
        results = [[],[]] # for monitoring progress
        name  = mir.split(':')[0]
        oType = classConvert[mir.split(':')[1]]
        if oTyp2use != oType: 
            continue
        
        # Get real data
        real_events = mirResultsPkl[name].matchEvents['A1_to_m8'][oType] 
        # Filter for AGAPs
        tmp_reals = []
        for event in real_events:
            for gene in event:
                if gene.startswith('AGAP'):
                    tmp_reals.append(gene)
                else:pass
        results[0].extend(tmp_reals)
        del(tmp_reals)
        
        #get ctrl data
        for i in range(len(mirResultsPkl[name].ctrlEvents['A1_to_m8'])):
            ctrl_events = mirResultsPkl[name].ctrlEvents['A1_to_m8'][i][oType]
            tmp_ctrls_i = []
            # filtr for agaps
            for event in ctrl_events:
                for gene in event:
                    if gene.startswith('AGAP'):
                        tmp_ctrls_i.append(gene)
                    else:pass
            results[1].append(tmp_ctrls_i)
            del(tmp_ctrls_i)
            
        
        # Add the single miR results to the master set
        masterDest[0].update(results[0])
        for i in range(len(results[1])):
            masterDest[1][i].update(results[1][i])
        None
        
    

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
                ha='center', va='bottom',rotation='vertical')
    



# --------- File Parsing ---------


anoXcsv = '/Users/biggus/Documents/James/Data/GeneOntology/Ag-Sep-2009-Web.csv'
A1_m8_file = '/Users/biggus/Documents/James/Data/Tu_miRNA/ResultsMiRNA_Targeting/2009_11_11/2009_11_11.AGAP.stats.100.medFDRmeth.A1-m8passed.txt'
pklPath_Ca = '/Users/biggus/Documents/James/Data/Tu_miRNA/ResultsMiRNA_Targeting/2009_11_11/2009_11_11.AGAP.seedMatches.100mvCtrls.storeEvents.pkl'
figTit     = ' title '
# -- AnoXcel --
print "\n\nprocessing AnoXcel..."
anoX_reader = csv.reader(open(anoXcsv, 'rU'))
anoXdata = {}
goInfo   = {}


domains = Bag({'bp':Bag({'start':137,
                         'end':142}),
               'mf':Bag({'start':127,
                         'end':132}),
               'cc':Bag({'start':132,
                         'end':137}),})

# Set Which Domain to look at
strt = domains.cc.start
stp  = domains.cc.end

# Set Event Class
useClass = 3


for row in anoX_reader:
    if not row[0].startswith('AGAP'):continue # dont want header
    
    if row[0][:-3] in anoXdata:
        anoXdata[row[0][:-3]].append([row[0]]+row[strt:stp])
        goInfo[row[strt:stp][-2]] = row[strt:stp] # Create a library of goTermInfo
    else:
        anoXdata[row[0][:-3]] = [[row[0]]+row[strt:stp]]
        goInfo[row[strt:stp][-2]] = row[strt:stp] # Create a library of goTermInfo


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
        real,stDv,med = fdrStats
        if med+stDv <= 0.5: # enforce FDR filter
            goFDRs.append((gTerm,real,med+stDv))
goFDRs.sort(key=lambda x: x[1])
goFDRs.reverse()
print "Printing GO FDRs..." 
for fdr in goFDRs:
    print "%s" % ('\t'.join([str(x) for x in fdr]))

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
ax.xaxis.set_minor_locator( mpl.ticker.AutoLocator() )

labelFDRs(rects,[x[2] for x in goFDRs])

pl.xticks(fontsize=9, rotation=90)
pl.show()