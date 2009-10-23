from Bio import SeqIO
import JamesDefs





#========================= User Defined Variables =========================

#  File paths: 
#  NOTE: genomeFileOne should get the seqs from the genome
#       that occurs FIRST in the orthologList file!!!

genomeFileOne = '/Users/biggus/Documents/James/Data/Tu_miRNA/Fastas/Cq_500afterCoding.newCoords.usuable.stpCdn.fas'
genomeFileTwo = '/Users/biggus/Documents/James/Data/Tu_miRNA/Fastas/Ag_500afterCoding.usuable.stpCdn.fas'
orthologList  = open('/Users/biggus/Documents/James/Data/OrthologDefs/Culex_Agam_1-to-1.nr.txt','r').readlines()

outFile       = '/Users/biggus/Documents/James/Data/Tu_miRNA/Fastas/OrthoPairFastas/CqAg_orthosForMDOS_500afterCoding.usuable.stpCdn.fas'

#==========================================================================




#  Instantiate the fasta rec lists with BioPython Seq using geneID field of discriptor as key to seq objects
genomeOneFastasDict = SeqIO.to_dict(SeqIO.parse(open(genomeFileOne, "rU"), 'fasta'),
                                    key_function = lambda rec : rec.description.split()[0])

genomeTwoFastasDict = SeqIO.to_dict(SeqIO.parse(open(genomeFileTwo, "rU"), 'fasta'),
                                    key_function = lambda rec : rec.description.split()[0])


#  Initiate resultList
resultList = []

#  Explode orthologList into list of lists
JamesDefs.explodeDelimitedList(orthologList, '\t')

#  Populate a list of GeneIDs in each genome's dict of boundary seqs
genomeOneGeneIDs = genomeOneFastasDict.keys()
genomeTwoGeneIDs = genomeTwoFastasDict.keys()

#  Loop through orthologList and call each fasta in orthoPair, format
#  the new comboFasta and append it to resultList
for orthoPair in orthologList:
    
    #  Test for orthoPair[0] in genomeOneFastasDict and same for orthoPair[1] in genomeTwoFastasDict
    orthoPair_0_warn = None
    orthoPair_1_warn = None
    if orthoPair[0] not in genomeOneGeneIDs:
        orthoPair_0_warn = 'Yes'
    if orthoPair[1] not in genomeTwoGeneIDs:
        orthoPair_1_warn = 'Yes'
    
    #  If either are not there, print warning
    if 'Yes' in [orthoPair_0_warn, orthoPair_1_warn]:
        
        notFoundList = []
        if orthoPair_0_warn == 'Yes':
            notFoundList.append(orthoPair[0])
        if orthoPair_1_warn == 'Yes':
            notFoundList.append(orthoPair[1])
            
        print 'There was a problem with orthoPair '+str(orthoPair)+'. The following were not found in expected lists: '+str(notFoundList)
        
    
    else:
        #  Format combo fasta such that the name is 'ortho1-ortho2' and seq is orthoSeq1+\n+orthoSeq2
        orthoComboFasta = '>'+orthoPair[0]+':'+orthoPair[1]+'\n'+genomeOneFastasDict[orthoPair[0]].seq.tostring()+'\n'+genomeTwoFastasDict[orthoPair[1]].seq.tostring()+'\n\n'
        resultList.append(orthoComboFasta)
        orthoComboFasta = ''


resultFile = open(outFile,'w')

resultFile.writelines(resultList)


print str(len(resultList))+' of '+str(len(orthologList))+' pairs where combined and included.\n'