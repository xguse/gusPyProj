

def getMultiGene(srtdTx,TxCol):
    """Removes recs whose gene occurs only once in super list."""
    srtdTx.sort(key=lambda x: x[TxCol])

    multiTx = []
    i=0
    lng = len(srtdTx)
    while i < lng:
        nxt = i+1
        pgrs = False
        while pgrs == False:
            iRec   = srtdTx[i]
            try:
                nxtRec = srtdTx[nxt]
            except:
                nxtRec = ['eof-eof-eof']*len(srtdTx[i])
            if iRec[TxCol][:-3] == nxtRec[TxCol][:-3]:
                nxt+=1
            else:
                pgrs = True
        if (nxt-i>1):
            multiTx.extend(srtdTx[i:nxt])
            i = nxt
        else:
            i = nxt
    return multiTx


def comboSortByGENE(pathA,pathB,colA,colB,oPath=None):
    """Given two tableFiles sort all rows by GeneName
    as deduced by TxName.  Must give col where TxName lives
    in each file. 
    If oPath: sorted table is also written out."""
    # +++++ read in data and add ID and Tx Cols +++++
    tblA = map(lambda l: [pathA.split('/')[-1],l.strip('\n').split('\t')[colA]] + l.strip('\n').split('\t'), open(pathA))
    tblB = map(lambda l: [pathB.split('/')[-1],l.strip('\n').split('\t')[colB]] + l.strip('\n').split('\t'), open(pathB))

    #tstA,tstB = tblA[0],tblB[0]

    combo = sorted(tblA+tblB,key=lambda x: x[1])
    if oPath:
        out = open(oPath,'w')
        for l in combo:
            out.write('%s\n' % ('\t'.join(l)))
        out.flush()
        out.close()
        return combo
    else:
        return combo


def hMMsplc_jDict2edgeBED(jDict,outPath):
    """Given a jDict with whole line, write info to minimal BEDfile
    using the edgeCoords as contigStrt/stop and no 'block' data."""
    oFile = open(outPath, 'w')
    # --- get list of all coverages ---
    for chrm in jDict:
        for jnc in jDict[chrm]:
            origLine = jDict[chrm][jnc][0].split('\t')
            oFile.write('%s\t%s\t%s\t%s\t%s\n' % (chrm, jnc[0], jnc[1], origLine[3], origLine[4]))

def hMMsplc_jDict2BED(jDict,outPath):
    """Given a jDict with whole line, write info to BEDfile."""
    oFile = open(outPath, 'w')
    # --- get list of all coverages ---
    for chrm in jDict:
        for jnc in jDict[chrm]:
            oFile.write('%s\n' % (jDict[chrm][jnc][0]))

