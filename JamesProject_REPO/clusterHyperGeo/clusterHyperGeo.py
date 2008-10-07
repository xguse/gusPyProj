from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
from Bio import SeqIO
from doug_hypergeometric import hyperGeoPvalue
import JamesDefs
import re
import string
from time import time

#--------- Script Specific Function Definitions ---------------------

def countMotifInAll(motifStr, seqDict):
    """
    Takes:	- motif string (RegEx string)
                - dict of seqObjs to search

    Does:	- compiles motifStr into regEx obj
                - resets each seq's motifPresent attrib to 0
                - uses regEx obj to loop over seqDict and update motifPresent attrib in each that seq contains motif				
    
    Returns:	- count of seqs containing motif and modifies seqObjs' motifPresent attrib in place
    """    
    
    import re
    
    
    #  Compile RegEx obj from motifStr
    motifRegEx = re.compile(motifStr, re.IGNORECASE)
    
    #  Loop over list and count motif in fwd and revComp oris for presence absense
    #  Sum total hits in totalHits
    
    totalHits = 0
    
    for record in seqDict:
        #  initiate motifPresent attrib to 0
        seqDict[record].motifPresent = 0
        
        if motifRegEx.search(seqDict[record].seq.tostring()) != None:
            # update motifPresent attrib to reflect found motif
            seqDict[record].motifPresent = 1
            totalHits += 1
    
        # check revComp only if above test did not find the motif
        if seqDict[record].motifPresent == 0:
            
            # check the revComp for motif and update record's motifPresent if found
            if motifRegEx.search(seqDict[record].seq.reverse_complement().tostring()) != None:
                # update motifPresent attrib to reflect found motif
                seqDict[record].motifPresent = 1
                totalHits += 1
    
    return totalHits

def convertMotifList(motifList):
    i = 0
    while i < len(motifList):
        motifList[i] = [motifList[i], JamesDefs.iupac2regex(motifList[i])]
        i += 1

#--------------------------------------------------------------------



#========================= User Defined Variables =========================
#  -Input Files-
clusterDefinitionList = map(string.strip, open('/Users/biggus/Documents/James/Data/Osvaldo/OM_MaleFemaleClusterDefs.txt', 'r'))
motifList             = map(string.strip, open('/Users/biggus/Documents/James/Data/2KB/2kb_Sequence/2kb_Combo/2Kb_AllMosquitoes/MosqMotifs/upstream_exclsv-conserved_mosquito-motifs_nr.txt', 'r'))
boundarySeqs          = '/Users/biggus/Documents/James/Data/2KB/2kb_Sequence/2kb_Anopheles/2KBupTSS_goodAffyAGAPsFastasOUT.masked.nr.fas'

#  -Output File-
outFile               = '/Users/biggus/Documents/James/Data/Osvaldo/OM_MaleFemale.mosqMotifs137.pVals.txt'
outFile = open(outFile,'w')

#==========================================================================

t1 = time()

#  Populate a Dict with Seq objs for Anopheles boundary seqs
#  What follows directly is a klugde to get my seqDict vals to have the IUPAC ambiguous alphabet
boundarySeqs = list(SeqIO.parse(open(boundarySeqs, "rU"), "fasta"))
for record in boundarySeqs :
    record.seq.alphabet = IUPACAmbiguousDNA

boundarySeqs = SeqIO.to_dict(boundarySeqs, key_function = lambda rec : rec.description.split()[0])

# convert iupac motifs to regexs and creat list of lists with each motif represented as ['IUPAC', 'REGEX'] 
convertMotifList(motifList)


#  group ClusterDefs by ClusterName
clusterDefinitionList = JamesDefs.groupByField(clusterDefinitionList, 0)

#  This will become a list of tab delim'd params for the hyperGeo func: 'Motif:ClusterID';'motifCountInAll';'len(all)';'motifCountInCluster';'numOfSeqsInCluster' 
hyperGeoParams_4_motifClusterPairs = []

m=0
for motif in motifList:
    m+=1
    print 'Motif '+str(m)
    #  Count how many seq in total list have motif in either orientation
    motifCountInAll = None
    motifCountInAll = countMotifInAll(motif[1], boundarySeqs)
 
    
    for cluster in clusterDefinitionList:

        #  count how many seq in Cluster have motifPresent set to 1
        motifCountInCluster = 0
        
        for seqID in cluster:
            if boundarySeqs.has_key(seqID[1]):
                motifCountInCluster += boundarySeqs[seqID[1]].motifPresent
            else:
                print seqID[1]+' not present in boundrySeq dictionary!'
            
        #  Format hyperGeoParams_4_motifClusterPairs record
        tab = '\t'
        # ------------------------------------->  no,n,ko,ki
        hyperGeoParams_4_motifClusterPairs.append(motif[0]+','+seqID[0]+tab+str(len(boundarySeqs))+tab+str(len(cluster))+tab+str(motifCountInAll)+tab+str(motifCountInCluster)) 

t2 = time()        
        
t3 = time()
c = 0
while c < len(hyperGeoParams_4_motifClusterPairs):
    #  explode string to list of params
    params = hyperGeoParams_4_motifClusterPairs[c].split('\t')
    
    #  calc hyperGeo pVal
    pVal = hyperGeoPvalue(int(params[1]),int(params[2]),int(params[3]),int(params[4]))
    print "%.6f" % (pVal)
    hyperGeoParams_4_motifClusterPairs[c] = hyperGeoParams_4_motifClusterPairs[c]+'\t'+str(pVal)+'\n'
    c+=1    
t4 = time()



outFile.writelines(hyperGeoParams_4_motifClusterPairs)




print 'The cluster counting took '+str((t2-t1)/60)+' minues.\nThe hyperGeo function took '+str((t4-t3)/60)+' minutes.'

