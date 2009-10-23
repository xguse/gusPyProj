from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
from Bio import SeqIO


#--------- Script Specific Function Definitions ---------------------


def countLENsAndNs(boundarySeqs):
    data = []
    
    
    for SeqID in boundarySeqs:
        #test = boundarySeqs[SeqID]
        data.append([boundarySeqs[SeqID].id, len(boundarySeqs[SeqID].seq), boundarySeqs[SeqID].seq.count('N')])
    print 'Counting done!'
    return data

def formatAndWrite(outFilePath, dataList):
    outFile = open(outFilePath, 'w')
    
    totals = [0,0]
    for boundaryRegion in dataList:
        totals[0] += boundaryRegion[1]
        totals[1] += boundaryRegion[2]
        
    outFile.write('#Total_bp:%i\tTotal_Ns:%i\tTotal_Non-Ns:%i\tTotal_Ns/bp:%.2f\n' % (totals[0],totals[1],totals[0]-totals[1],(float(totals[1])/totals[0])*100))
    
    print 'Totals done!'
    print '#Total_bp:%i\tTotal_Ns:%i\tTotal_Non-Ns:%i\tTotal_Ns/bp:%.4f\n' % (totals[0],totals[1],totals[0]-totals[1],(float(totals[1])/totals[0])*100)
    
    outFile.write('GeneID\tbp\tNs\tNonNs\tNs/bp\n')
        
    for i in dataList:
        geneID  = i[0]
        bp      = i[1]
        Ns      = i[2]
        NonNs   = bp-Ns
        percent = (float(Ns)/bp)*100
        
        outFile.write('%s\t%i\t%i\t%i\t%.2f\n' % (geneID,bp,Ns,NonNs,percent))
    
    print 'All done!'
    outFile.close()
        

#--------------------------------------------------------------------



#========================= User Defined Variables =========================
boundarySeqs = '/Users/biggus/Documents/James/Data/2KB/2kb_Sequence/2kb_Aedes/aedes2KBupStreamTSS.masked.fas'

outFilePath  = '/Users/biggus/Documents/James/Data/2KB/2kb_Sequence/2kb_Aedes/aedes2KBupStreamTSS.N_analysis.txt'

#==========================================================================




#  Populate a Dict with Seq objs for boundary seqs
#  What follows directly is a klugde to get my seqDict vals to have the IUPAC ambiguous alphabet
boundarySeqs = list(SeqIO.parse(open(boundarySeqs, "rU"), "fasta"))
for record in boundarySeqs :
    record.seq.alphabet = IUPACAmbiguousDNA

boundarySeqs = SeqIO.to_dict(boundarySeqs, key_function = lambda rec : rec.description.split()[0])

#test = boundarySeqs[boundarySeqs.keys()[0]]


#  will hold a list of seq ID, length, and number of Ns per seq as secondary lists 
#  Ex: [[Seq1,389,45],[Seq2,1790,100]]
dataList = countLENsAndNs(boundarySeqs)

#  format and write data to file
formatAndWrite(outFilePath, dataList)





