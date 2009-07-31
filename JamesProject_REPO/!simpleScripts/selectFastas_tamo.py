from TAMO.seq import Fasta
import string





#========================= User Defined Variables =========================

#  Path to original file
##originalFastaDict = open('/Users/biggus/Documents/James/Data/2KB/2kb_Sequence/2kb_Anopheles/2KBupTSS_goodAffyAGAPsFastasOUT.masked.nr.fas', 'rU')

##desiredFastaList  = open('/Users/biggus/Documents/James/Data/ClusterDefs/TC-46.txt', 'rU').readlines()

##outFile           = '/Users/biggus/Documents/James/Data/2KB/2kb_Sequence/2kb_Combo/2Kb_AllMosquitoes/MosqMotifs/MotifGroupPWMs/AllGroups/ModuleData/TC-46.fas'

originalFastaDict = '/Users/biggus/Documents/James/Data/2KB/2kb_Sequence/2kb_Aedes/aedes2KBupStreamTSS.softMasked.geneStrand.fas'

desiredFastaList  = open('/Users/biggus/Documents/James/Collaborations/Campbell/data/CCupAt4Days.genes.txt', 'rU').readlines()

outFile           = '/Users/biggus/Documents/James/Collaborations/Campbell/data/CCupAt4Days.UNmasked.fas'

hardMask          = None

#==========================================================================

#  Strip newlines from fasta ID list
desiredFastaList = map(string.strip, desiredFastaList)


#  Instantiate the fasta rec lists
originalFastaDict = Fasta.load(originalFastaDict)

#  New dict to catch copied seqObjs
desiredFastaDict = {}

for rec in desiredFastaList:
    if originalFastaDict.has_key(rec):    
        desiredFastaDict[rec] = originalFastaDict[rec]
    else:
        print rec+' not found in source fasta list!'
    

# Hard Mask if requested
if hardMask:
    for x in desiredFastaDict:
        desiredFastaDict[x] = desiredFastaDict[x].replace('a','N')
        desiredFastaDict[x] = desiredFastaDict[x].replace('c','N')
        desiredFastaDict[x] = desiredFastaDict[x].replace('g','N')
        desiredFastaDict[x] = desiredFastaDict[x].replace('t','N')
else:
    for x in desiredFastaDict:
        desiredFastaDict[x] = desiredFastaDict[x].upper() # make sure all letters are uppercase for downstream compatibility

#  Write selected recs to outFile

Fasta.write(desiredFastaDict, outFile)

print "Done."