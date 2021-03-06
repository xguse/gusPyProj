from Bio.Blast import NCBIXML
import cPickle

# === <defs> === #

    
def writeOutMatchData(blastObjs,outFileHandle,contigStrech):
    """
    Takes blastObjs, an outFileHandle and the chosen shortest contiguous length
    of alignments to report. Writes out each kept probes list of overlaps with its
    matching probes.
    """
    overlapGroups = {}
    recs = blastObjs.keys()
    recs.sort()
    for rec in recs:
        blastRec = blastObjs[rec]
        olpGroup = set()
        for a in range(len(blastRec.alignments)):
            for h in range(len(blastRec.alignments[a].hsps)):
                if blastRec.alignments[a].hsps[h].frame[0] == -1:  # only look at matches in SAME orientation bc RevComp of mRNA is used to hybridize.
                    continue
                longestMatchSeq = len(max(blastRec.alignments[a].hsps[h].match.split()))
                alnPcntID       = float(blastRec.alignments[a].hsps[h].identities)/blastRec.alignments[a].hsps[h].align_length*100
                if longestMatchSeq >= contigStrech:
                    olpGroup.add(blastRec.alignments[a].hit_id)
                    outFileHandle.write('%s\t%s\t%s\t%.1f\t%s\t%s\t%s\t%s\t%s\t%s\n' % \
                        (blastRec.query_id, 
                         blastRec.alignments[a].hit_id,
                         blastRec.alignments[a].hsps[h].align_length,
                         alnPcntID,
                         longestMatchSeq,
                         blastRec.alignments[a].hsps[h].query_start,
                         blastRec.alignments[a].hsps[h].query_end,
                         blastRec.alignments[a].hsps[h].sbjct_start,
                         blastRec.alignments[a].hsps[h].sbjct_end,
                         blastRec.alignments[a].length))
        overlapGroups[blastRec.query_id] = list(olpGroup)
     
    # Write out each probe's match list at a threshold of contigStrech
    outFileHandle.write("\n%s\n\nBelow are the probes that match more than themselves with at least %s bp of contiguous identities:\n" % (' - '*20,contigStrech))
    for rec in recs:
        assert rec in overlapGroups[rec], "ERROR: overlapGroups[%s] does not contain itself!"
        if len(overlapGroups[rec]) > 1:  # only write out sets if they match more than theirself.
            overlapGroups[rec].sort()
            outFileHandle.write('%s:\t%s\n' % (rec,'\t'.join(overlapGroups[rec])))


            
# === </defs> === #



# === <main> === #
xmlFile = '/Users/biggus/Documents/James/Data/MicroArrayPrep/Aedes/Probes/ProbeBlastResults/Aedes.agilent.probes.nr.blastn_low.AedesTxs.xml'
usePkl  = '/Users/biggus/Documents/James/Data/MicroArrayPrep/Aedes/Probes/ProbeBlastResults/Aedes.agilent.probes.nr.blastn_low.AedesTxs.pkl'
bObjPkl = '/Users/biggus/Documents/James/Data/MicroArrayPrep/Aedes/Probes/ProbeBlastResults/Aedes.agilent.probes.nr.blastn_low.AedesTxs.pkl'
outFile = '/Users/biggus/Documents/James/Data/MicroArrayPrep/Aedes/Probes/ProbeBlastResults/Aedes.agilent.probes.nr.blastn_low.AedesTxs.15stretch.txt'

#xmlFile = '/Users/biggus/Documents/James/Data/MicroArrayPrep/Aedes/Probes/ProbeBlastResults/Aedes.agilent.probes.nr.test.blastn_low.self.xml'
#outFile = '/Users/biggus/Documents/James/Data/MicroArrayPrep/Aedes/Probes/ProbeBlastResults/Aedes.agilent.probes.nr.test.blastn_low.self.redunProbes.txt'

contigStretch = 15

if not usePkl:
    blastObjs      = {}
    
    for rec in NCBIXML.parse(open(xmlFile)):
        blastObjs[rec.query_id] = rec
        # Pickle blastObjs if file path is given and we are not using a pkl
        if bObjPkl:
            print 'Pickle to be written to %s' % (bObjPkl)
            cPickle.dump(blastObjs,open(bObjPkl,'w'))
else:
    blastObjs      = cPickle.load(open(usePkl,'rU'))    

# Write out a list of each probes alignments that has a contiguous streatch of identities
# at least as long as contigStretch.  Then write out a nr list of each probe's matching Txs.
outFile = open(outFile, 'w')



# write out header info
outFile.write('All Alignments >= %s considered.\n' % (contigStretch))
outFile.write('Query_Probe\tMatching_Tscript\tLength_Of_Alignment\tPcntID_Of_Alignment\tLongest_Contiguous_Run\tProbe_Start\tProbe_End\tTscript_Start\tTscript_End\tTscript_Length\n')
    
writeOutMatchData(blastObjs,outFile,contigStretch)

    
print 'Done.'
# === </main> === #


toDo = \
"""_done_ -> finish writing writeOutMatchData()"""