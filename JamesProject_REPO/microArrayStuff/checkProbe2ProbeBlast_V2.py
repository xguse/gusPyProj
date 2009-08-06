from Bio.Blast import NCBIXML

# === <defs> === #
def retainOrNot2(Bio_Blast_Record_Blast,contigStrech):# full name of BioPython obj for clairity
    """Determines whether blast record should be kept or passed over.
    Returns record if it is to be retained and a Set obj of identicle
    probe_ids. tupple:(blastRec, setOfProbeIds)"""

    blastRec = Bio_Blast_Record_Blast # renamed for ease of use
    
    # Is there more than one alignment object?
    if len(blastRec.alignments) < 2:
        return None
    
    # Collect ids of all alignments with an hsp of 100% ident
    # along given stretch of the query seq in same orientation.
    idents = []
    for a in range(len(blastRec.alignments)):
        for h in range(len(blastRec.alignments[a].hsps)):
            if blastRec.alignments[a].hsps[h].frame[0] == -1:  # make sure the probes are in same orientation
                continue
            
            longestMatchSeq = len(max(blastRec.alignments[a].hsps[h].match.split(' ')))
            if longestMatchSeq >= contigStrech:
                idents.append(blastRec.alignments[a].title)
                
    idents = frozenset(idents)
    
    if len(idents) > 1:
        return (blastRec,idents)
    else:
        return None
                
    
def writeOutMatchData(group,keptBlastObjs,outFileHandle,contigStrech):
    """Takes group of matching probes, keptBlastObjs, and an outFileHandle.
    Writes out each kept probes list of overlaps with its matching probes."""
    
    for i in range(len(group)):
        blastRec = keptBlastObjs[group[i].strip()]
        
        for a in range(len(blastRec.alignments)):
            for h in range(len(blastRec.alignments[a].hsps)):
                if blastRec.alignments[a].hsps[h].frame[0] == -1:
                    continue
                if blastRec.alignments[a].hit_id == blastRec.query_id:
                    continue
                longestMatchSeq = len(max(blastRec.alignments[a].hsps[h].match.split(' ')))
                if longestMatchSeq >= contigStrech:
                    outFileHandle.write('%s\t%s\t%s\t%.3f\t%s\t%.3f\t%s\n' % \
                        (blastRec.query_id, \
                         blastRec.alignments[a].hit_id, \
                         blastRec.query_letters,
                         float(blastRec.alignments[a].hsps[h].identities)/blastRec.query_letters,\
                         longestMatchSeq,\
                         float(getLongestHspStretch(blastRec.alignments[a]))/blastRec.query_letters,\
                         getLongestHspStretch(blastRec.alignments[a]),))


def getLongestHspStretch(BioBlastAlignObj):
    longestStrech = 0
    
    for hsp in BioBlastAlignObj.hsps:
        if len(max(hsp.match.split())) > longestStrech:
            longestStrech = len(max(hsp.match.split()))
        return longestStrech
            
# === </defs> === #



# === <main> === #
xmlFile = '/Users/biggus/Documents/James/Data/MicroArrayPrep/Aedes/Probes/ProbeBlastResults/Aedes.agilent.probes.blastn_low.AedesTxs.xml'
outFile = '/Users/biggus/Documents/James/Data/MicroArrayPrep/Aedes/Probes/ProbeBlastResults/Aedes.agilent.probes.blastn_low.AedesTxs.15stretch.test.txt'

#xmlFile = '/Users/biggus/Documents/James/Data/MicroArrayPrep/Aedes/Probes/ProbeBlastResults/Aedes.agilent.probes.nr.test.blastn_low.self.xml'
#outFile = '/Users/biggus/Documents/James/Data/MicroArrayPrep/Aedes/Probes/ProbeBlastResults/Aedes.agilent.probes.nr.test.blastn_low.self.redunProbes.txt'

length = 15

# Sort out which probes are identicle and create a dict of
# their blastRecs and a set of each group of redundant probes.
keptBlastObjs   = {}
matchingProbes = set([])
for rec in NCBIXML.parse(open(xmlFile)):
    recTuple = retainOrNot2(rec,length)
    if recTuple != None:
        keptBlastObjs[recTuple[0].query_id] = recTuple[0]
        matchingProbes.add(recTuple[1])

# Write out the list of identicle probe groups followed a list
# of their alignments.
outFile = open(outFile, 'w')

'''may want to write out a lot of source file info right up front'''
        
matchingProbes = list(matchingProbes)
for i in range(len(matchingProbes)):
    matchingProbes[i] = list(matchingProbes[i])
    matchingProbes[i].sort()
matchingProbes.sort()

outFile.write('All Alignments >= %s considered.\n%s matching groups.\n' % (length,len(matchingProbes)))

for group in matchingProbes:
    outFile.write('%s\n' % ('\t'.join(group)))

outFile.write('\n%s\n\n' % (' - '*20)) 
# write out header info
outFile.write('All Alignments >= %s considered.\n' % (length))
outFile.write('Query_Probe\tMatching_Probe\tQuery_Length\tFraction_of_ID(HSP)\tLongest_Stretch_of_ID(HSP)\tLongest_Fraction_of_ID(MatProbe)\tLongest_Stretch_of_ID(MatProbe)\n')
    
for group in matchingProbes:
    writeOutMatchData(group,keptBlastObjs,outFile,length)

    
print 'Done.'
# === </main> === #


toDo = \
"""_done_ -> finish writing writeOutMatchData()"""