import cPickle
from gusPyCode.defs.JamesDefs import DotDict
from gusPyCode.defs.statsDefs import benjHochFDR
from gusPyCode.defs.bioDefs import goEnrichment

def goClassEnrichment(geneCluster,geneClusterName,goClassDict,popSize=None,FDRthresh=0.05):
    """
    goClassEnrichment(geneCluster,goClassDict,popSize,FDRthresh=0.05):
    geneCluster     = set(genesGroupedBySomeQuality)
    geneClusterName = UniqeID
    goClassDict     = dict(keys=each GO term In bp, cc, or mf class, vals=set(genesAttchd2GOterm)
    popSize         = int(numberOfGenesConsideredAsPopulation)
    
    Returns List for each GeneSet GOterm Combo:
    [GOterm,GOTermSize,geneSetName,geneClstSize,pVal,BHpVal,FDRthresh,numMatchingGenes,matchingGenes].
    """
    rawPs = []
    for goTerm in goClassDict:
        rawPs.append([goTerm,goEnrichment(geneCluster,goClassDict[goTerm],popSize=popSize)])
    adjPs = benjHochFDR(rawPs)
    
    for i in range(len(adjPs)):
        gTrm      = adjPs[i][0]
        p         = adjPs[i][1]
        bhP       = adjPs[i][2]
        mtchGenes = goClassDict[gTrm].intersection(geneCluster)
        adjPs[i]  = [gTrm,
                     str(len(goClassDict[gTrm])),
                     geneClusterName,
                     str(len(geneCluster)),
                     '%.5g' %(p),
                     '%.5g' %(bhP),
                     str(FDRthresh),
                     str(len(mtchGenes)),
                     str(sorted(list(mtchGenes)))]
    
    
    return adjPs

def filtAGAP(listOfTups):
    rList =[]
    for i in range(len(listOfTups)):
        for j in range(len(listOfTups[i])):
            if listOfTups[i][j].startswith('AGAP'):
                rList.append(listOfTups[i][j])
    return set(rList)
    
    
if __name__ == '__main__':
    import psyco
    psyco.profile()
    
    statsFile = map(lambda l: l.strip('\n').split(' : '), open('/Users/biggus/Documents/James/Data/Tu_miRNA/ResultsMiRNA_Targeting/2009_10_26run/2009_10_26.AGAP.100bothCtrls.stats.events.medFDRmeth.txt','rU').readlines())
    goPkl = DotDict(cPickle.load(open('/Users/biggus/Documents/James/Data/GeneOntology/goStats.pkl','r')))
    

    for i in range(10,13):
        if statsFile[i][1].startswith('all'):
            print '|'.join(statsFile[i])
            bpGO = goClassEnrichment(filtAGAP(eval(statsFile[i][-1])),'%s_%s'%(statsFile[i][0],statsFile[i][1]),goPkl.bp2Genes,popSize=12457)
            for each in bpGO:
                if float(each[4]) < float(each[6]):
                    print ' >%s'%('\t'.join(each))
        else:
            print '|'.join(statsFile[i])
            bpGO = goClassEnrichment(filtAGAP(statsFile[i][-1]),'%s_%s_%s'%(statsFile[i][0],statsFile[i][1],statsFile[i][2]),goPkl.bp2Genes,popSize=12457)
            for each in bpGO:
                if float(each[4]) < float(each[6]):
                    print ' >%s'%('\t'.join(each))
        
    