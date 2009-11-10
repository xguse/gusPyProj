from Bio import SeqIO
import string
from gusPyCode.defs import JamesDefs




#========================= User Defined Variables =========================

#  Path to original file
originalFastaDict = open('/Users/biggus/Documents/James/Data/2KB/2kb_Sequence/2kb_Anopheles/2KBupTSS_goodAffyAGAPsFastasOUT.masked.nr.fas', 'rU')

desiredFastaList  = '/Users/biggus/Documents/James/Data/ClusterDefs/TC-Clusters.txt'

outDir            = '/Users/biggus/Documents/James/Data/ClusterDefs/TC-Fastas/'

#==========================================================================

desiredFastaList = map(lambda line : line.strip(), open(desiredFastaList, 'rU').readlines())

# Parse clusterDefs into list of clusters
listOfClusterDefs = JamesDefs.groupByField(desiredFastaList,0)




#  Instantiate the fasta rec lists with BioPython Seq using geneID field of discriptor as key to seq objects
originalFastaDict = SeqIO.to_dict(SeqIO.parse(originalFastaDict, 'fasta'),
                                    key_function = lambda rec : rec.description.split()[0])

for cluster in listOfClusterDefs:
    print "Working on Cluster: %s" % (cluster[0][0])
    #  New dict to catch copied seqObjs
    desiredFastaObjList = []
    
    for rec in cluster:
        if originalFastaDict.has_key(rec[1]):    
            desiredFastaObjList.append(originalFastaDict[rec[1]])
        else:
            print rec[1]+' not found in source fasta list!'
        
    #  Write selected recs to outFile
    
    outFile = "%s%s.fas" % (outDir, cluster[0][0])
    outFile = open(outFile, 'w')
    SeqIO.write(desiredFastaObjList, outFile, 'fasta')
    outFile.close()

print "Done."