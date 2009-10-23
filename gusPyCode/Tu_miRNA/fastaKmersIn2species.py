from TAMO.seq import Fasta
import pickle
# User Vars #
orthoTabMatch = '/Users/biggus/Documents/James/Data/Tu_miRNA/SeedCountOutPut/counts/AaCq_matureSeedTabs.combined.txt'
orthoTabCtrl  = '/Users/biggus/Documents/James/Data/Tu_miRNA/SeedCountOutPut/ctrls/AaCq_matureSeedTabs.combined.ctrl.txt'
##kmerFasta     = '/Users/biggus/Documents/James/Data/Tu_miRNA/MDOSoutPut/AaCq_7mers_seedHitsFromMDOS.motifzSgte3.ctrl.txt'
oFile         = '/Users/biggus/Documents/James/Data/Tu_miRNA/SeedCountOutPut/targets/AaCq_seedTargets.2genomes.miRB.txt'
pFile         = '/Users/biggus/Documents/James/Data/Tu_miRNA/SeedCountOutPut/pickles/AaCq_seedTargets.2genomes.miRB.countPickle.txt'
#-----------#

sourceFiles = '# matchFile = %s\n# ctrlFile = %s\n# outFile = %s\n# Data Pickle = %s\n\n' % \
            (orthoTabMatch.split('/')[-1], \
             orthoTabCtrl.split('/')[-1],\
             oFile.split('/')[-1],\
             pFile.split('/')[-1])

orthoTabMatch = open(orthoTabMatch, 'rU').readlines()
#print orthoTabMatch[0],'\n\n',orthoTabMatch[1]
headInfoMatch  = [orthoTabMatch.pop(0),orthoTabMatch.pop(0)][1].strip().split('\t')[1:]
orthoTabMatch = ''.join(orthoTabMatch).split('--')

orthoTabCtrl = open(orthoTabCtrl, 'rU').readlines()
headInfoCtrl  = [orthoTabCtrl.pop(0),orthoTabCtrl.pop(0)][1].strip().split('\t')[1:]
orthoTabCtrl = ''.join(orthoTabCtrl).split('--')

# check head infos
hdInfoMtCk = []
hdInfoClCk = []
for i in headInfoMatch:
    hdInfoMtCk.append(i.split('_')[0])
for i in headInfoCtrl:
    hdInfoClCk.append(i.split('_')[0])
assert hdInfoMtCk == hdInfoClCk, "Header info does not seem to match between the orthoTab files."


    
    

# convert "string\nstring\n" into [list,list] after all is done
# __Match group__
for i in range(0,len(orthoTabMatch)):
    # convert to [string,string]
    ##x = orthoTabMatch[i]
    orthoTabMatch[i] = orthoTabMatch[i].strip().split('\n')
    ##x = orthoTabMatch[i]
    for j in range(0, len(orthoTabMatch[i])):
        # convert to [list,list]
        orthoTabMatch[i][j] = orthoTabMatch[i][j].split('\t')
# __Ctrl group__        
for i in range(0,len(orthoTabCtrl)):
    # convert to [string,string]
    orthoTabCtrl[i] = orthoTabCtrl[i].strip().split('\n')
    
    for j in range(0, len(orthoTabCtrl[i])):
        # convert to [list,list]
        orthoTabCtrl[i][j] = orthoTabCtrl[i][j].split('\t')

        
resultsMtch = {}
for each in headInfoMatch:
    resultsMtch[each]=[]
resultsCtrl = {}
for each in headInfoCtrl:
    resultsCtrl[each]=[]

# Determine which miRNA seed if any has hits in both orthologs
# Count "hit" as either hit in pos1 OR pos2
# __Match group__
for orthoPair in orthoTabMatch:
    for i in range(1,len(orthoPair[0])):
        if orthoPair[0][i] == '1' and orthoPair[1][i] == '1':
            resultsMtch[headInfoMatch[i-1]].append('%s:%s' % (orthoPair[0][0],orthoPair[1][0]))
# __Ctrl group__
for orthoPair in orthoTabCtrl:
    for i in range(1,len(orthoPair[0])):
        if orthoPair[0][i] == '1' and orthoPair[1][i] == '1':
            resultsCtrl[headInfoCtrl[i-1]].append('%s:%s' % (orthoPair[0][0],orthoPair[1][0]))

rKeysMtch = resultsMtch.keys()
rKeysMtch.sort()
rKeysCtrl = resultsCtrl.keys()
rKeysCtrl.sort()

seedCounts = []
miRNA_MtchCtrl = {'!doc':['matchHits','CtrlHits']}
for i in range(len(rKeysMtch)):
    seedCounts.append('%s:%s\n' % (rKeysMtch[i], len(resultsMtch[rKeysMtch[i]])))
    seedCounts.append('%s:%s\n' % (rKeysCtrl[i], len(resultsCtrl[rKeysCtrl[i]])))
    seedCounts.append('-- --\n')
    miRNA_MtchCtrl['_'.join([rKeysMtch[i].split('_')[0], rKeysMtch[i].split('_')[-1]])] = [len(resultsMtch[rKeysMtch[i]]),len(resultsCtrl[rKeysCtrl[i]])] 
    #print '# %s:%s\n' % (rKeysMtch[i], len(resultsMtch[rKeysMtch[i]]))
    #print '# %s:%s\n' % (rKeysCtrl[i], len(resultsCtrl[rKeysCtrl[i]]))
#seedCounts.append('# ** ** ** ** **\n\n')
#seedCounts.append('''# What follows are the orthologous genes with seed
## hits (miRNA[pos1-7]OR[pos2-8]) to miRNAs indicated.\n\n''')

#
##oList = ['# Counts of Seeds provided:\n']
##oList.extend(seedCounts)
##for each in rKeys:
    ##eList = ['%s\n' % (each)]
    ##for e in results[each]:
        ##eList.append('%s\n'%(e))
    ##eList.append('--\n\n')
    ##oList.extend(eList)
    
oFile = open(oFile, 'w')
pFile = open(pFile, 'w')
oFile.write(sourceFiles)
oFile.writelines(seedCounts)
pickle.dump(miRNA_MtchCtrl,pFile)


print 'done'