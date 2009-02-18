from TAMO import MotifMetrics,MotifTools
from string import maketrans

seqFile       = '/Users/biggus/Documents/James/Data/2KB/2kb_Sequence/2kb_Aedes/aedes2KBupStreamTSS.softMasked.geneStrand.fas'
seqsOfIntrest = '/Users/biggus/Documents/James/Writings_Talks/Grants/09_Feb/POI_Aa.txt'
motifFile     = '/Users/biggus/Documents/James/Data/2KB/2kb_Sequence/2kb_Combo/2Kb_AllMosquitoes/MosqMotifs/upstream_exclsv-conserved_mosquito-motifs_nr.txt'

# print out files used
print 'seqFile: %s' % (seqFile)
print 'seqsOfIntrest: %s' % (seqsOfIntrest)
print 'motifFile: %s' % (motifFile)

seqDB         = MotifMetrics.Fasta.file2dict(seqFile)
seqsOfIntrest = map(lambda line: line.strip(), open(seqsOfIntrest,'rU').readlines())
motifList     = map(lambda line: line.strip(), open(motifFile,'rU').readlines())

# initiate translation table for converting softMask to Hard mask
trans = maketrans('acgt','NNNN')





# convert list of motifs into list of lists(each primary index == [motif,TAMO_motif_obj])
for i in range(0, len(motifList)):
    motifList[i] = [motifList[i], MotifTools.Motif_from_text(motifList[i])]
  
# scan each seqOfIntrest for each motif and report how many times each occurs 

for i in seqsOfIntrest:
    tSeq = seqDB[i].translate(trans)  # translate will convert any softmasked seqs to HardMasked (scan does not ignor lower)
    results = ['----------\nSeq:%s len:%s Ns:%s\nMotif\tInstances_Found\tStarting_Positions_of_Instances'% (i,len(tSeq), tSeq.count('N'))]
    for j in motifList:
        hits = j[1].scan(tSeq)  
        if hits[0]:
            startPos = []
            for p in hits[1]:
                startPos.append(str(p[0]))
            results.append('%s\t%s\t%s' % (j[0],len(hits[0]),','.join(startPos)))
    
    # print out resuts for in in seqsOfIntrest
    for line in results:
        print line
    print ''