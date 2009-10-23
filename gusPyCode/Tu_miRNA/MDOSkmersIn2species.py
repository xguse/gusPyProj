

# User Vars #
orthoTabs = '/Users/biggus/Documents/James/Data/Tu_miRNA/SeedCountOutPut/ctrls/AaCq_seedTabs.combined.ctrl.txt'
zScrsFile = '/Users/biggus/Documents/James/Data/Tu_miRNA/MDOSoutPut/AaCq_7mers_seedHitsFromMDOS.motifzSgte3.ctrl.txt'
oFile     = '/Users/biggus/Documents/James/Data/Tu_miRNA/SeedCountOutPut/targets/AaCq_seedTargets.2genomes.ctrl.txt'
#------------

orthoTabs = open(orthoTabs, 'rU').readlines()
headInfo  = [orthoTabs.pop(0),orthoTabs.pop(0)][1].strip().split('\t')[1:]
orthoTabs = ''.join(orthoTabs).split('--')
 
zScrsFile = map(lambda line: line.strip(), open(zScrsFile,'rU').readlines())
zScrsGTE3 = []
for z in zScrsFile:
    if z.startswith('>'):
        zScrsGTE3.append(z[1:])
    
    

# convert "string\nstring\n" into [list,list] after all is done
for i in range(0,len(orthoTabs)):
    # convert to [string,string]
    orthoTabs[i] = orthoTabs[i].strip().split('\n')
    
    for j in range(0, len(orthoTabs[i])):
        # convert to [list,list]
        orthoTabs[i][j] = orthoTabs[i][j].split('\t')
#o = orthoTabs[0]
        
results = {}
for each in headInfo:
    results[each]=[]

# Determine which miRNA seed if any has hits in both orthologs
for orthoPair in orthoTabs:
    for i in range(1,len(orthoPair[0])):
        if orthoPair[0][i] == '1' and orthoPair[1][i] == '1':
            results[headInfo[i-1]].append('%s:%s' % (orthoPair[0][0],orthoPair[1][0]))

rKeys = results.keys()
rKeys.sort()
seedCounts = []
for each in rKeys:
    seedCounts.append('# %s:%s\n' % (each, len(results[each])))
seedCounts.append('# -- --\n\n')
seedCounts.append('''# What follows are the orthologous genes "targeted" by the seed (miRNA[pos2-8]) indicated.
# Targets of seeds that were assigned conservation z-scores of >= 3 by MDOS are reported 
# here with stars.\n\n''')

oList = ['# Counts of Seeds provided:\n']
oList.extend(seedCounts)
for each in rKeys:
    if each in zScrsGTE3:
        eList = ['%s *\n' % (each)]
    else:
        eList = ['%s\n' % (each)]
    for e in results[each]:
        eList.append('%s\n'%(e))
    eList.append('--\n\n')
    oList.extend(eList)
    
oFile = open(oFile, 'w')
oFile.writelines(oList)

print 'done'