

def countMotifInAll(motifStr, seqDict):
    """
    Takes:	- motif string (RegEx string)
                - Dict of Bio.Seq objects with SeqID as Keys

    Does:	- processes BioSeqDict to produce a bdrySeqList
                - passes motifStr and bdrySeqList to countMotif()
                - receives count of seqs w/ motif
                            
    Returns:	- count of seqs in total list of seqs containing motif
    """
    
    secList = []
    
    for k, v in seqDict:
        secList.append(v.seq.tostring())
        
        
        
# genomeTwoFastasDict[orthoPair[1]].seq.tostring()