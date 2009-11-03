import numpy
from TAMO import seq
from TAMO.Clustering import MotifCompare
from TAMO import MotifTools
from TAMO.MotifTools import Motif
from TAMO.MD.Meme import Meme
from TAMO.MD.AlignAce import AlignAce
from gusPyCode.defs.JamesDefs import randFromList_noReplace
from gusPyCode.defs import xpermutations
import random

def alignPairWithOffSet(motif1, motif2, offset):
    """
    Takes two TAMO motifs and the offset for motif2 to motif1 (can be negative).
    Returns a string representing the text alignment dilimited by '\n'.
    """

    line1 = None
    line2 = None

    
    if offset < 0:
        line1 = '%s%s' % ('-' * abs(offset),motif1.oneletter)
        line2 = motif2.oneletter
    elif offset > 0:
        line1 = motif1.oneletter
        line2 = '%s%s' % ('-' * abs(offset),motif2.oneletter)
    else:
        line1 = motif1.oneletter
        line2 = motif2.oneletter
    
    # right pad if needed   
    if len(line1) != len(line2):
        if len(line1) == min([len(line1),len(line2)]):
            line1 = '%s%s' % (line1,'-' * abs(len(line1)-len(line2)))
        else:
            line2 = '%s%s' % (line2,'-' * abs(len(line1)-len(line2)))
        
    return '%s\n%s' % (line1,line2)


def alignAndCombineMotifs(motifs, weights):
    # zip motifs and weights
    simMotifs = zip(motifs, weights)
    # sort by weights
    simMotifs.sort(key=lambda x: abs(x[1]))
    simMotifs.reverse()
    
    aligned = alignSimilarMotifs([x[0] for x in simMotifs], minoverlap=4)
    #print '--'
    #for each in aligned: print each.oneletter
    #print '\n'
    
    comboMotif = MotifTools.sum(aligned,[-x[1] for x in simMotifs])
    return comboMotif


def findBestPairAlignments(listOfMotifObjs, minoverlap=6, verbose=None):
    """
    Takes: list of TAMO motif objects.  Finds best pairwise alignments among list members, trying both
    orientations. Motifs in list are numbered by original index in results. Returns: 2D list of 
    results for each combination of motifs with the matrix coords corresponding to motif 
    index in original list (exp: dist of motif0 and motif4 == 2dList[0][4]; BUT 2dList[4][0] == None).
    Always put lower index first or you will get 'None'.  Same index twice also gives 'None'
    
    Each value at the 2D coords contains a tuple: (alignOri,distScore,alignment,offset).
    alignOri = 1 == both motifs in original ori.  alignOri = -1 == motif with higher index
    was revComped to get best score.
    
    verbose == True prints the scores, orientations and alignments for each motif pair.
    """
    # rename listOfMotifObjs for brevity
    motifs = listOfMotifObjs
    
    # Initialize empty return-matrix
    rMat = []
    for i in range(len(motifs)):
        rMat.append([None]*len(motifs))
        
    # Create list of non-redundant index combos for comparing
    toCompare = [x for x in xpermutations.xuniqueCombinations(range(len(motifs)),2)]
    
    for i in range(len(toCompare)):
        alignOri  = None
        distScore = None
        alignment = None
        offset    = None
        
        minDiffOri = getMinDiffOri(motifs[toCompare[i][0]],motifs[toCompare[i][1]],minoverlap=minoverlap, getOffset=1)
        
        # If pos ori, then motif obj returned will be ref to original motifs[toCompare[i][1]]
        # else: newly constructed revComp is returned
        if motifs[toCompare[i][1]] is minDiffOri[0]: alignOri = 1
        else: alignOri = -1
        
        distScore = minDiffOri[1]
        alignment = alignPairWithOffSet(motifs[toCompare[i][0]], minDiffOri[0], minDiffOri[2])
        offset = minDiffOri[2]
        
        # Assign tuple to  matrix coords:
        rMat[toCompare[i][0]][toCompare[i][1]] = (alignOri, distScore, alignment, offset)

    # Write out the results if verbose
    if verbose:
        oString = '#MotifPair\tAlignOri\tAlignScore\tAlignment\n'
        for pair in toCompare:
            tmp = '%s:%s\t%s\t%.3G\t%s' \
                % (str(pair[0])+'_'+motifs[pair[0]].oneletter,
                   str(pair[1])+'_'+motifs[pair[1]].oneletter,
                   rMat[pair[0]][pair[1]][0],
                   rMat[pair[0]][pair[1]][1],
                   rMat[pair[0]][pair[1]][2])
            # Futz with formating to allow alignments to match when pasted in an exclFile
            tmp    = tmp.split('\n')
            spc    = ' '*3
            add    = spc.join(['\t']*tmp[0].count('\t'))
            tmp[1] = '%s%s%s\n' % (spc,add,tmp[1])
            
            print '\n'.join(tmp)
            oString += '\n'.join(tmp)
            
        
    
    return rMat
    
    


# ----------------------
def getKmersWithOneMisMtch(motif1, motifListWithMetrics):
    """
    Takes a TAMO motif and a list of lists: [TAMOmotif, weightMetric].
    Returns listOfLists in same form, if a revComp in motifListWithMetrics
    matches better, _IT_ is returned instead of the original motif.
    
    Restricts motif1 from collecting any motifs that are not of leng +/- 1 of itself.
    """
    resultList = []
    
    for mWithMetric in motifListWithMetrics:
        if len(motif1)-1 <= len(mWithMetric[0]) <= len(motif1)+1:
            # Determine what distanceResult == one misMatch.
            # Use length of shortest motif for misMatch Calc
            maxAlnLen   = min(len(motif1), len(mWithMetric[0]))
            kMer        = Motif('%s' %('A'*maxAlnLen))
            kMer1mis    = Motif('%s%s' %('A'*(maxAlnLen-1),'T'))  # whether at end or in middle the 'T gives same align score'
            oneMisMatch = MotifCompare.minshortestoverhangdiff(kMer,kMer1mis)
            bestOri = getMinDiffOri(motif1,mWithMetric[0])
            if bestOri[1] <= oneMisMatch:
                resultList.append([bestOri[0],mWithMetric[1]]) # keep [bestOriTAMOmotif, weightMetric]
            
    return resultList


# ----------------------
def alignSimilarMotifs(similarMotifs, minoverlap=6):
    """
    Takes list of similar motifs. List should be sorted by weight or pvalue since this def
    basically aligns each motif to the top motif.  Its kind of a progressive pairwise alignment.
    """
    alignedMotifs = []
    matrix = findBestPairAlignments(similarMotifs, minoverlap=minoverlap, verbose=None)
    
    # Find longest neg-offset to motif0
    lNegOff = 0
    for i in range(1,len(similarMotifs)):
        if matrix[0][i][3] < lNegOff:
            lNegOff = matrix[0][i][3]
    # Adjust left padding of motif0 for longest negOffset
    alignedMotifs.append(similarMotifs[0][lNegOff,similarMotifs[0].width])
    None
    
    # Adjust orientation and left padding for each motif
    for i in range(1,len(similarMotifs)):
        if matrix[0][i][0] < 0:                                   # If neg offset
            rcMotif = similarMotifs[i].revcomp()                  #   revComp motif_i 
            lPad = lNegOff-matrix[0][i][3]                        #   remember -> neg - neg = closer to 0
            alignedMotifs.append(rcMotif[lPad,rcMotif.width])
        elif matrix[0][i][0] > 0:                                 # If pos offset
            lPad = lNegOff-matrix[0][i][3]                        #   remember -> neg - pos = farther from 0
            alignedMotifs.append(similarMotifs[i][lPad,similarMotifs[i].width])
        elif matrix[0][i][0] == 0:                                # If no offset
            alignedMotifs.append(similarMotifs[i][lNegOff,similarMotifs[i].width])
    
    # Add right padding to match longest motif from above
    # Find Longest motif
    lMotifLen = 0
    for mtf in alignedMotifs:
        if mtf.width > lMotifLen:
            lMotifLen = mtf.width
    for i in range(len(alignedMotifs)):
        if alignedMotifs[i].width < lMotifLen:
            alignedMotifs[i] = alignedMotifs[i][0,lMotifLen]
        
        
    
    return alignedMotifs


# ----------------------
def trimPaddedMotif(paddedMotif):  # doesnt seem needed since motifTools.sum() seems to do this automatically
    """
    Takes a padded motif and removes blank positions from left and right, until it
    sees a non-blank.  Returns trimmed motif obj.
    """
    
    # Count Blanks
    lBlanks     = 0
    rBlanks     = 0
    posInfoBits = paddedMotif.bits[:] # real copy bc i will be playing with orders and dont want to cause trouble
    for pos in posInfoBits:
        if pos == 0:
            lBlanks+=1
        else:
            break
    posInfoBits.reverse()   
    for pos in posInfoBits:
        if pos == 0:
            rBlanks+=1
        else:
            break
        
    # Remove Blanks
    motif = paddedMotif[lBlanks:-rBlanks+1]
    
    return motif

# ----------------------
def getMinDiffOri(motif1,motif2,minoverlap=6,getOffset=False):     ##originally had this at end of func def. dunno why->   , N=1, keepLen=0):
    """
    Takes two TAMO motifs.  Calculates TAMO.Clustering.MotifCompare.minshortestoverhangdiff
    for motif1 against motif2 and the rvcmp of motif2.  Returns a tuple containing the TAMO
    motif obj of motif2 that produced the least distance result and the distance result.
    
    (motif2, distResult) -OR- (motif2_rc, distResult)
    if getOffset:
    (motif2, distResult, offset) -OR- (motif2_rc, distResult, offset)
    """
    
    dist = MotifCompare.minshortestoverhangdiff(motif1,motif2,minoverlap=minoverlap,want_DistAndOff=1)

        
    if getOffset:    
        if dist[2]:
            return (motif2.revcomp(), dist[0],dist[1])
        else:
            return (motif2, dist[0],dist[1])
    else:
        if dist[2]:
            return (motif2.revcomp(), dist[0])
        else:
            return (motif2, dist[0])
    
# ----------------------
def genRandClusters(geneList,totalSeqs, N=1, keepLen=0):
    """
    Takes geneList and returns N non-overlapping lists of randomly chosen geneNames of
    len(geneList) from totalSeqs. If keepLen != 0, length of individual sequences are
    kept bt +/- 5% of original sequences.
    <totalSeqs can be filePath OR sequence Dictionary>
    """
    
    # NOTE TO SELF: may want to look at using FakeFasta for this in the future.
    # For Now:
    # Open totalSeqs and load into dict if not a dict.
    if type(totalSeqs) != type({}):
        totalSeqs = seq.Fasta.file2dict(totalSeqs)
        
    
    clustSize = len(geneList)
    
    # build new lists of geneNames
    if keepLen:
        randClusts = []
        for i in range(N):
            tempClus = []
            for gene in geneList:
                geneLen = len(totalSeqs[gene])
                fivePercent = geneLen*0.025
                allGenes = totalSeqs.keys()
                
                while allGenes:
                    tempGene = randFromList_noReplace(allGenes)
                    if len(totalSeqs[tempGene]) < geneLen+fivePercent:
                        if len(totalSeqs[tempGene]) > geneLen-fivePercent:
                            if tempGene not in tempClus:
                                inRandClus = False
                                for l in randClusts:
                                    if tempGene in randClusts:
                                        inRandClus = True
                                        continue
                                if inRandClus == False:
                                    tempClus.append(tempGene)
                                    break
            # Warn if you got all the way through allGenes
            # without fililing up tempCluster to geneList
            # length.
            assert len(tempClus) == len(geneList), 'tempClus_%s did not reach len(geneList)!'
            randClusts.append(tempClus)
        
                                    
    else:
        # Place code that doesnt care about seq len here
        print "This part has NOT been Validated!!\n"*5
        randClusts = []
        for i in range(N):
            tempClus = []
            for gene in geneList:
                allGenes = totalSeqs.keys()
                
                while allGenes:
                    tempGene = randFromList_noReplace(allGenes)
                    if tempGene not in tempClus:
                        inRandClus = False
                        for l in randClusts:
                            if tempGene in randClusts:
                                inRandClus = True
                                continue
                        if inRandClus == False:
                            tempClus.append(tempGene)
                            break
            # Warn if you got all the way through allGenes
            # without fililing up tempCluster to geneList
            # length.
            assert len(tempClus) == len(geneList), 'tempClus_%s did not reach len(geneList)!'
            randClusts.append(tempClus)
        
    return randClusts


# ----------------------
def transfacLike2tamoMotif(filePath,beta=0.01,bg={'A':.25,'C':.25,'G':.25,'T':.25}): 
    """
    Takes file path of a single matrix right now.  Returns a tuple with tuple[0] = name
    and tuple[1] = TAMO motif object initiated with the source matrix.
    
    beta and bg defaults are set to those in TAMO.MotifTools.Motif_from_counts()
    """

    
    tfMatrix = map(lambda l: l.strip(), open(filePath).readlines())
    name     = ''
    pwm      = []
    
    # this version only expects one matrix per file
    numOfRecs = 0
    for l in tfMatrix:
        if l.startswith('DE'):
            numOfRecs+=1
    assert numOfRecs == 1, \
           'Error: transfacLike2tamoMotif() expects only one "DE" line. %s found.' % (numOfRecs)
    
    
    for l in tfMatrix:
        if l.startswith('DE'):
            # take the first tab separated field that is not 'DE' as the name
            name = l.split('\t')[1]
            break  ## only ONE name please
    
    # read string version of PWM into pwm
    for l in tfMatrix:
        try:
            type(int(l.split('\t')[0]))
        except:
            continue
        
        lineType = type(int(l.split('\t')[0]))
        if lineType == type(0):
            pwm.append(l.split('\t')[1:5])
    
    # number-ify pwm
    for i in range(len(pwm)):
        # create dict of Nuc values
        pwm[i] = {'A':float(pwm[i][0]),'C':float(pwm[i][1]),'G':float(pwm[i][2]),'T':float(pwm[i][3])}

    tamoMotif = MotifTools.Motif_from_counts(pwm,beta=beta,bg=bg)
    return (name,tamoMotif)

# ----------------------
def shuffleSeqDict(seqDict):
    """Takes dict of seqs and returns a dict with identicle keys but with each corresponding
DNA seq randomly shuffled with location of Ns preserved."""
    
    assert type(seqDict) is dict, 'shuffleSeqDict() takes exactly one dictionary.'
    shuffledDict = {}
    
    # Will remove Ns, shuffle remaining nucleotides, then insert back into the
    # native N environment.
    for k in seqDict:
        seqList     = list(seqDict[k])
        seqNucsList = list(seqDict[k].replace('N',''))       
        random.shuffle(seqNucsList)
        random.shuffle(seqNucsList)
        nuc = 0
        for i in range(0,len(seqList)):
            if seqList[i] != 'N':
                seqList[i] = seqNucsList[nuc]
                nuc += 1
        shuffledDict[k] = ''.join(seqList)
        
    return shuffledDict

# ----------------------
def loadMotifsFromOutFile(fileName,MD_outType=None,):
    validOptions = ['Meme', 'AlignAce', 'MDscan', 'Weeder','list']
    assert MD_outType in validOptions, 'MD_outType=%s; valid options are:\n%s' \
           % (MD_outType, '\n'.join(validOptions))
    if MD_outType == validOptions[0]:
        # read meme lines into TAMO MEME class and parse for motifs
        fileNameLines = map(lambda l: l.strip('\n'), open(fileName, 'rU').readlines())
        
        # Initiate MD class, read in lines, parse lines, and extract motifs
        memeObj = Meme()
        memeObj.lines = fileNameLines
        memeObj._parse()
        return memeObj.motifs
    
    # fill out other MD apps
    elif MD_outType == validOptions[1]:
        # read meme lines into TAMO MEME class and parse for motifs
        fileNameLines = map(lambda l: l.strip('\n'), open(fileName, 'rU').readlines())
        
        # Initiate MD class, read in lines, parse lines, and extract motifs
        aceObj = AlignAce()
        aceObj.lines = fileNameLines
        aceObj._parse()
        return aceObj.motifs
    
    elif MD_outType == validOptions[4]:
        # read lines into list
        fileNameLines = map(lambda l: l.strip('\n'), open(fileName, 'rU').readlines())
        
        # convert text motifs into TAMO motif objs
        for i in range(len(fileNameLines)):
            fileNameLines[i] = MotifTools.Motif_from_text(fileNameLines[i])
        return fileNameLines


# ----------------------
def seqSubSet(seqNames,sourceFastaPath):
    """Takes list of seq names and a fasta file.
    Returns tuple -> (dictOfSubSet, listOfMissingSeqNames)"""
    
    # create dict of source Fasta File
    allSeqs = seq.Fasta.load(sourceFastaPath)
    
    subSet = {}
    missingNames =[]
    for i in seqNames:
        if allSeqs[i]:
            subSet[i] = allSeqs[i]
        else:
            missingNames.append(i)
    
    return (subSet, missingNames)


changeLog = """2009-04-16 -- Established ChangeLog
2009-04-16 -- added randClusterGen
2009-06-05 -- added getMinDiffOri
"""
