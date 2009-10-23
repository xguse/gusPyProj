import cPickle
from gusPyCode.defs.JamesDefs import initList
from gusPyCode.defs import mathDefs
from gusPyCode.defs import miRNA_targeting_V05 as miRT


print '\n\n'

outFile = '/Users/biggus/Documents/James/Data/Tu_miRNA/ResultsMiRNA_Targeting/testPickle.5randMiRs.100bothCtrls.stats.events.medFDRmeth.prgData.txt'#'/Users/biggus/Documents/James/Data/Tu_miRNA/ResultsMiRNA_Targeting/2009_08_31run.100permCtrls.stats.events.indSeeds.medFDRmeth.txt'
outFile = open(outFile, 'w')



pklPath_Ca = '/Users/biggus/Documents/James/Data/Tu_miRNA/ResultsMiRNA_Targeting/testPickle.5randMiRs.100prmCtrls.strEvents.pkl'#'/Users/biggus/Documents/James/Data/Tu_miRNA/ResultsMiRNA_Targeting/2009_08_31.seedMatches.100permCtrls.pkl' #2009_09_29.seedMatches.100psCtrls.pkl'
data_Ca    = cPickle.load(open(pklPath_Ca,'rU'))

pklPath_Cb = '/Users/biggus/Documents/James/Data/Tu_miRNA/ResultsMiRNA_Targeting/testPickle.5randMiRs.100psCtrls.strEvents.pkl'#'/Users/biggus/Documents/James/Data/Tu_miRNA/ResultsMiRNA_Targeting/2009_08_31.seedMatches.100permCtrls.pkl' #2009_09_29.seedMatches.100psCtrls.pkl'
data_Cb    = cPickle.load(open(pklPath_Cb,'rU'))

consFdrThreshold = 0.25
stdvsAboveMed    = 2

# calc the various centers and spreads
items = ['aga-miR-1890',
         'aga-miR-9c',
         'aga-miR-13b',
         'aga-miR-133',
         'aga-miR-11',
         'aga-miR-988',
         'aga-miR-278',
         'aga-miR-965',]

#miR_matches = {}
#for i in items:
    #miR_matches[i] = data['miR_matches'][i]

miR_matches_Ca = data_Ca['miR_matches']
miR_matches_Cb = data_Cb['miR_matches']

assert sorted(miR_matches_Ca.keys()) == sorted(miR_matches_Cb),\
       'It seems your Ca and Cb miRNAs do not match.'

# purge redun data
print 'Purging objs...'
for obj in miR_matches_Ca:
    miR_matches_Ca[obj].purgeData()
    miR_matches_Cb[obj].purgeData()


print 'Getting Gene Targets...'


##def writeTargetsOldMeth(miRobj,oFile):
    ##""" WARNING: Parameters set above (consFdrThreshold,stdvsAboveMed) do not pertain to this def
                 ##You must set them IN THE DEF."""
    ##miRHits = miR_matches[miR].reportGeneTargets(zScore=3,FDR=0.15)
        
    ##totals   = (None,[],[],[])  # (None,[geneNumFor1,meanFDRFor1,genesIn1], [sameFor2], [sameFor3]

    
    ### Get totaled data
    ##for orthoType in range(1,4):
        ##print 'mir: %s oType: %s' % (miR_matches[miR].name,orthoType)
        ##totGenes = set()
        ##FDRs     = []
        ##for seedType in miRHits:
            ##print seedType
            ##try:
                ##if miRHits[seedType][orthoType]:
                    ##totGenes.update(miRHits[seedType][orthoType][0])
                    ##FDRs.append([len(miRHits[seedType][orthoType][0]),miRHits[seedType][orthoType][2]]) # [numOfGenes,FDR]
            ##except:
                ##exit('an error occured some where in "try:" on line48.')
            
                
        ##if len(totGenes) > 0:
            ##comboFDR = sum([float(x[0])*x[1] for x in FDRs])/sum([x[0] for x in FDRs])
            ##totals[orthoType].extend([len(totGenes),comboFDR,sorted(list(totGenes))])
        ##else:
            ##totals[orthoType].extend([0,None ,['None']])
    
    ### write out totals data
    ##for i in range(1,len(totals)):
        
        ##outFile.write('%s : allPassedSeedsFor_%s : %s : %s / %s\n' \
                      ##%(miR_matches[miR].name,
                        ##i,
                        ##totals[i][0],
                        ##totals[i][1],
                        ##','.join(sorted([str(x) for x in totals[i][2]]))))
    ### write out passed Seed data
    ##for seedType in sorted(miRHits):
        ##for i in range(1,len(miRHits[seedType])):
            ##if miRHits[seedType][i]:
                ##outFile.write('%s : %s : orthoType_%s : %s : %.2f_c : %.2f_m / %s\n'\
                              ##%(miR_matches[miR].name,
                                ##seedType,
                                ##i,
                                ##len(miRHits[seedType][i][0]),
                                ##miRHits[seedType][i][2],
                                ##miRHits[seedType][i][1],
                                ##','.join(sorted([str(x) for x in miRHits[seedType][i][0]]))))
    ##outFile.write('-- --\n')
    ##outFile.flush()

def writeTargetsFdrMedMeth(miRobj_Ca,miRobj_Cb,oFile):
    #print 'Processing %s...' % (miRobj_Ca.name)
    miRHits_Ca = miRobj_Ca.reportGeneTargetsFdrMedMeth(stdvLimit=stdvsAboveMed,
                                                       consFdrThresh=consFdrThreshold,
                                                       divide=0)
    totReal = [None,set(),set(),set()]
    totCtrl = initList(len(miRobj_Ca.ctrlEvents[miRobj_Ca.ctrlEvents.keys()[0]]),[None, set(), set(),set()])
    
    # Calulate combined FDR for miRNA using Ctrl_b data from seedTypes that passed
    # the reportGeneTargetsFdrMedMeth() Ctrl_a screen.
    #   >> Gather and combine data from passed seedTypes:
    for oType in range(1,4):
        for sType in miRT._seedModels:
            if miRHits_Ca[sType][oType]:
                trLen_0 = len(totReal[oType])
                totReal[oType].update(miRHits_Ca[sType][oType][0])
                rUpdtLen = len(miRHits_Ca[sType][oType][0])
                trLen_1 = len(totReal[oType])
                None
                for i in range(len(totCtrl)):
                    tciLen_0 = len(totCtrl[i][oType])
                    totCtrl[i][oType].update(miRobj_Cb.ctrlEvents[sType][i][oType])
                    cUpdtLen = len(miRobj_Cb.ctrlEvents[sType][i][oType])
                    tciLen_1 = len(totCtrl[i][oType])
                    None
            
    #   >> Calculate separate FDRs for each Ctrl_b group:
    totalsData = [None,None,None,None]
    for oType in range(1,4):
        if totReal[oType] == set():
                continue
        tempFDRs = []
        for i in range(len(totCtrl)):
            ctrlVal = len(totCtrl[i][oType])
            realVal = len(totReal[oType])
            
            if float(ctrlVal)/realVal >= 1:
                tempFDRs.append(1.0)
            else:
                tempFDRs.append(float(ctrlVal)/realVal)
        tLen = len(tempFDRs)
        oFDRstdv,oFDRmed = mathDefs.stdDv(tempFDRs,'median')
        cons_oFDR        = oFDRmed + (stdvsAboveMed*oFDRstdv)
        totalsData[oType] = [totReal[oType],oFDRmed,cons_oFDR]
    
    # Write out Totals data:
    print miRobj_Cb.name
    outFile.write('-- %s --\n' % (miRobj_Cb.name)) 
    for i in range(1,len(totalsData)):
        if totalsData[i]:
            outFile.write('%s : allPassedSeedsFor_%s : %s : %s : %s  Seqs=%s\n' \
                          %(miRobj_Cb.name,
                            i,
                            len(totalsData[i][0]),
                            totalsData[i][2],
                            totalsData[i][1],
                            ','.join(sorted([str(x) for x in totalsData[i][0]]))))

                    
        
    # write out passed Seed data
    for seedType in sorted(miRHits_Ca):
        for i in range(1,len(miRHits_Ca[seedType])):
            if miRHits_Ca[seedType][i]:
                outFile.write('%s : %s : orthoType_%s : %s : %.2f : %.4f Seqs=%s\n'\
                              %(miRobj_Ca.name,
                                seedType,
                                i,
                                len(miRHits_Ca[seedType][i][0]),
                                miRHits_Ca[seedType][i][2],
                                miRHits_Ca[seedType][i][1],
                                ','.join(sorted([str(x) for x in miRHits_Ca[seedType][i][0]]))))
    outFile.flush()


assert sorted(miR_matches_Ca.keys()) == sorted(miR_matches_Cb.keys()),\
       "Your Ca and Cb files do not seem to contain the same miRNAs."
for miR in miR_matches_Ca.keys():
    #writeTargetsOldMeth(miR_matches[miR],outFile)
    writeTargetsFdrMedMeth(miR_matches_Ca[miR],miR_matches_Cb[miR],outFile)
    
        
        

    
print 'Im Done.'