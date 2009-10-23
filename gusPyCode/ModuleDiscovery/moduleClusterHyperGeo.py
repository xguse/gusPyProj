from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
from Bio import SeqIO
from doug_hypergeometric import hyperGeoPvalue
import JamesDefs
from defs_moduleClusterHyperGeo import *
import re
import string
from time import time




#========================= User Defined Variables =========================
#  -Input Files-
clusterDefinitionList = map(string.strip, open('/Users/biggus/Documents/MBGB/Rotations/James/Data/ClusterDefs/052708_3xTukeyHSD_tissueClusters.txt', 'r'))
motifList             = map(string.strip, open('/Users/biggus/Documents/MBGB/Rotations/James/Data/testData4Mapper/MGandUpAt24_70percent.motifs', 'r'))
boundarySeqs          = '/Users/biggus/Documents/MBGB/Rotations/James/Data/2KB/2kb_Sequence/2kb_Anopheles/2kb_anophelesProcessed/2KB_anophelesUpstream/2KBupTSS_goodAffyAGAPsFastasOUT.masked.nr.fas'

mostMotifsInSet       = 3  # how many motifs you want in largest module?


#  -Output File-
outFile               = '/Users/biggus/Documents/MBGB/Rotations/James/Data/testData4Mapper/MGandUpAt24_70percent.pVals.txt'
outFile = open(outFile,'w')

#==========================================================================

t1 = time()

#  Populate a Dict with Seq objs for Anopheles boundary seqs
#  What follows directly is a klugde to get my seqDict vals to have the IUPAC ambiguous alphabet
boundarySeqs = list(SeqIO.parse(open(boundarySeqs, "rU"), "fasta"))
for record in boundarySeqs :
    record.seq.alphabet = IUPACAmbiguousDNA

boundarySeqs = SeqIO.to_dict(boundarySeqs, key_function = lambda rec : rec.description.split()[0])

##  I think that I will compile each motif at time of use this time to facilitate passing the combo as a list
# convert iupac motifs to regexs and create list of lists with each motif represented as ['IUPAC', 'REGEX'] 
#convertMotifList(motifList)

#  Generate a list of all combinations of spplied motifs that are 
listOfSearchModules = createSearchModules(motifList, 3)


#  group ClusterDefs by ClusterName
clusterDefinitionList = JamesDefs.groupByField(clusterDefinitionList, 0)

#  This will become a list of tab delim'd params for the hyperGeo func: 'Module:ClusterID';'moduleCountInAll';'len(all)';'moduleCountInCluster';'numOfSeqsInCluster' 
hyperGeoParams_4_moduleClusterPairs = []

m=0
for module in listOfSearchModules:
    m+=1
    
    moduleString = str(module) #  for later label use
    
    
    modTimeStrt = time()

    #  Count how many seq in total list have module in either orientation
    moduleCountInAll = None
    moduleCountInAll = countModuleInAll(module, boundarySeqs)
 
    
    for cluster in clusterDefinitionList:

        #  count how many seq in Cluster have modulePresent set to 1
        moduleCountInCluster = 0
        
        for seqID in cluster:
            if boundarySeqs.has_key(seqID[1]):
                moduleCountInCluster += boundarySeqs[seqID[1]].modulePresent
            else:
                print seqID[1]+' not present in boundrySeq dictionary!'
            
        #  Format hyperGeoParams_4_moduleClusterPairs record
        tab = '\t'
        # ------------------------------------->  no,n,ko,ki
        hyperGeoParams_4_moduleClusterPairs.append(moduleString+','+seqID[0]+tab+str(len(boundarySeqs))+tab+str(len(cluster))+tab+str(moduleCountInAll)+tab+str(moduleCountInCluster)) 
        len_hyperGeoParams_4_moduleClusterPairs= len(hyperGeoParams_4_moduleClusterPairs)
    modTimeStop = time()
    print 'Module %s %s took %s seconds.' %  (str(m), moduleString, str(modTimeStop-modTimeStrt))
    
t2 = time()        
        
t3 = time()
c = 0
while c < len(hyperGeoParams_4_moduleClusterPairs):
    #  explode string to list of params
    params = hyperGeoParams_4_moduleClusterPairs[c].split('\t')
    
    #  calc hyperGeo pVal
    pVal = hyperGeoPvalue(int(params[1]),int(params[2]),int(params[3]),int(params[4]))
    print pVal
    hyperGeoParams_4_moduleClusterPairs[c] = hyperGeoParams_4_moduleClusterPairs[c]+'\t'+str(pVal)+'\n'
    c+=1    
t4 = time()



outFile.writelines(hyperGeoParams_4_moduleClusterPairs)




print 'The cluster counting took '+str((t2-t1)/60)+' minues.\nThe hyperGeo function took '+str((t4-t3)/60)+' minutes.'

