from Bio import SeqIO
import string





#========================= User Defined Variables =========================

#  Path to original file
##originalFastaDict = open('/Users/biggus/Documents/James/Data/2KB/2kb_Sequence/2kb_Anopheles/2KBupTSS_goodAffyAGAPsFastasOUT.masked.nr.fas', 'rU')

##desiredFastaList  = open('/Users/biggus/Documents/James/Data/ClusterDefs/TC-46.txt', 'rU').readlines()

##outFile           = '/Users/biggus/Documents/James/Data/2KB/2kb_Sequence/2kb_Combo/2Kb_AllMosquitoes/MosqMotifs/MotifGroupPWMs/AllGroups/ModuleData/TC-46.fas'

originalFastaDict = open('/Users/biggus/Documents/James/Data/2KB/2kb_Sequence/2kb_Aedes/aedes2KBupStreamTSS.softMasked.geneStrand.fas', 'rU')

desiredFastaList  = open('/Users/biggus/Documents/James/Collaborations/Campbell/data/CCupAt4Days.genes.txt', 'rU').readlines()

outFile           = '/Users/biggus/Documents/James/Collaborations/Campbell/data/CCupAt4Days.masked.fas'

hardMask          = False 

#==========================================================================

#  Strip newlines from fasta ID list
desiredFastaList = map(string.strip, desiredFastaList)


#  Instantiate the fasta rec lists with BioPython Seq using geneID field of discriptor as key to seq objects
originalFastaDict = SeqIO.to_dict(SeqIO.parse(originalFastaDict, 'fasta'),
                                    key_function = lambda rec : rec.description.split()[0])

#  New dict to catch copied seqObjs
desiredFastaObjList = []

for rec in desiredFastaList:
    if originalFastaDict.has_key(rec):    
        desiredFastaObjList.append(originalFastaDict[rec])
    else:
        print rec+' not found in source fasta list!'
    

# Hard Mask if requested
if hardMask:
    for x in desiredFastaObjList:
        desiredFastaObjList[x] = desiredFastaObjList[x].replace('a','N')
        desiredFastaObjList[x] = desiredFastaObjList[x].replace('c','N')
        desiredFastaObjList[x] = desiredFastaObjList[x].replace('g','N')
        desiredFastaObjList[x] = desiredFastaObjList[x].replace('t','N')

#  Write selected recs to outFile

outFile = open(outFile, 'w')
SeqIO.write(desiredFastaObjList, outFile, 'fasta')
outFile.close()

print "Done."