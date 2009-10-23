from Bio.Blast import NCBIXML

# === <defs> === #
def retainOrNot(Bio_Blast_Record_Blast):# full name of BioPython obj for clairity
    """Determines whether blast record shoudl be kept or passed over.
    Returns record if it is to be retained and a Set obj of identicle
    probe_ids. tupple:(blastRec, setOfProbeIds)"""

    blastRec = Bio_Blast_Record_Blast # rename for ease of use
    
    # Is there more than one alignment object?
    if len(blastRec.alignments) < 2:
        return None
    
    # Collect ids of all alignments with an hsp of 100% ident
    # along 100% of the query seq in same orientation.
    idents = []
    for a in range(len(blastRec.alignments)):
        for h in range(len(blastRec.alignments[a].hsps)):
            if blastRec.alignments[a].hsps[h].frame[0] == -1:
                continue
            
            longestMatchSeq = len(max(blastRec.alignments[a].hsps[h].match.split(' ')))
            if blastRec.query_letters == longestMatchSeq:
                idents.append(blastRec.alignments[a].title)
                
    idents = frozenset(idents)
    
    if len(idents) > 1:
        return (blastRec,idents)
    else:
        return None
                
def writeAlignments(identProbes,keptBlastObjs,outFileHandle):
    """Takes group of identicle probes and writes their alignments to outfile."""
    blastRec = keptBlastObjs[identProbes[0].strip()]
    
    outFileHandle.write('Alignments for group(%s):\n' % ('; '.join(group)))

    for a in range(len(blastRec.alignments)):
        for h in range(len(blastRec.alignments[a].hsps)):
            if blastRec.alignments[a].hsps[h].frame[0] == -1:
                continue
            if blastRec.alignments[a].hit_id == blastRec.query_id:
                continue
            
            longestMatchSeq = len(max(blastRec.alignments[a].hsps[h].match.split(' ')))
            if blastRec.query_letters == longestMatchSeq:
                #outFileHandle.write('Alignments for group(%s):\n' % ('; '.join(group)))
                outFileHandle.write('\t%s:%s\n' % (blastRec.query_id, blastRec.alignments[a].hsps[h].query))
                outFileHandle.write('\t%s %s\n' % (' '*len(blastRec.query_id), blastRec.alignments[a].hsps[h].match))
                outFileHandle.write('\t%s:%s\n\n' % (blastRec.alignments[a].hit_id, blastRec.alignments[a].hsps[h].sbjct))
    
    outFileHandle.write('%s\n\n' % ('='*90))
# === </defs> === #



# === <main> === #
xmlFile = '/Users/biggus/Documents/James/Data/MicroArrayPrep/Aedes/Probes/ProbeBlastResults/Aedes.agilent.probes.nr.blastn_low.self.xml'
outFile = '/Users/biggus/Documents/James/Data/MicroArrayPrep/Aedes/Probes/ProbeBlastResults/Aedes.agilent.probes.nr.blastn_low.self.redunProbes.txt'


# Sort out which probes are identicle and create a dict of
# their blastRecs and a set of each group of redundant probes.
keptBlastObjs   = {}
identicleProbes = set([])
for rec in NCBIXML.parse(open(xmlFile)):
    recTuple = retainOrNot(rec)
    if recTuple != None:
        keptBlastObjs[recTuple[0].query_id] = recTuple[0]
        identicleProbes.add(recTuple[1])

# Write out the list of identicle probe groups followed a list
# of their alignments.
outFile = open(outFile, 'w')

'''will want to write out a lot of source file info right up front'''
        
identicleProbes = list(identicleProbes)
for i in range(len(identicleProbes)):
    identicleProbes[i] = list(identicleProbes[i])
    identicleProbes[i].sort()
identicleProbes.sort()

for group in identicleProbes:
    outFile.write('%s\n' % ('\t'.join(group)))

outFile.write('\n%s\n\n' % (' - '*20)) 
    
for group in identicleProbes:
    writeAlignments(group,keptBlastObjs,outFile)

    
print 'Done.'
# === </main> === #