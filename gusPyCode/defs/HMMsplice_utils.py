#!/usr/local/bin/python
#
# initialUtils.py
#
# Author: Michelle Dimon, May 2009
#

from types import StringTypes
import random
import string
import hmmErrors

##################################################################################################
#  Utilities to work with input files (fasta, fastq, etc)
##################################################################################################

def convertQualityStr(strVersion, offset=33):
    """Converts a quality string to a list of integers.
    an offset of 33 is Phred style and an offet of 64 is Solexa style"""
    return map( lambda x: ord(x)-offset, strVersion )

def convertToQualityStr(intList, offset=33):
    """Converts a list of integers to a quality string."""
    return "".join( map(lambda x: chr(x+offset), intList) )

def qualityToProbability(qual, isSolexa=True):
    """Given a score, returns the probability associated with that score.
    This is the probability that the base call is an error, so a very low number
    (like 10 ^ -5) is good quality.
    Calculations are done differently for solexa and sanger, default is solexa"""

    if isSolexa:
        return 1 / (1 + 10 ** (qual/10.0))
    else:
        return 10 ** (-qual / 10.0)

def reverseComplement(seq):
    """return the reverse complement of a sequence.
    Original code from Kael and Dale."""
    seq=seq.upper()
    # complement
    compl = complement(seq)
    # reverse
    return compl[::-1]

def complement(seq,transl=None):
    """Return complement of seq.
    Original code from Kael and Dale.
    """
    transl = string.maketrans('aAcCgGtTnNxX-\t\n ','tTgGcCaAnNxX-\t\n ')
    compl = seq.translate(transl)
    return compl

def randomlySelectFromFile(inputFile, outputFile, numToSelect, numLines):
    """Randomly selects 'numToSelect' lines from inputFile and puts them
    in outputFile.
    """

    # pick 'numToSelect' numbers from range(0, numLines)
    linesToUse = random.sample(xrange(numLines), numToSelect)

    linesDict = {}
    for lineNum in linesToUse:
        linesDict[lineNum] = 0

    # read through inputFile and if line number is in random set, write to output file
    count = 0
    out = open(outputFile, "w")
    for line in open(inputFile):     
        if linesDict.has_key(count):
            out.write(line)
        count += 1

    out.close()

def randomlySelectFromFastqFile(inputFastq, outputFile, numToSelect, numRecords, startPos=None, stopPos=None, useFasta=False):
    """Randomly selects 'numToSelect' fastq records from inputFile and puts them
    in outputFile.
    """

    # pick 'numToSelect' numbers from range(0, numLines)
    linesToUse = random.sample(xrange(numRecords), numToSelect)

    linesDict = {}
    for lineNum in linesToUse:
        linesDict[lineNum] = 0

    # read through inputFile and if line number is in random set, write to output file
    count = 0
    out = open(outputFile, "w")
    if useFasta:
        for (title, seq) in FastaIterator(inputFastq):
            if linesDict.has_key(count):
                if startPos != None:
                    seq = seq[startPos:]
                if stopPos != None:
                    seq = seq[:stopPos]
                out.write(">%s\n%s\n" % (title, seq))
            count += 1
    else:
        for (title, seq, qual) in FastqIterator(inputFastq):
            if linesDict.has_key(count):
                if startPos!= None:
                    seq = seq[startPos:]
                if stopPos != None:
                    seq = seq[:stopPos]
                out.write("@%s\n%s\n+%s\n%s\n" % (title, seq, title, qual))
            count += 1


    out.close()

def convertFastqToFasta(inputFastq, outputFasta):
    """Converts a fastq file into a fasta file."""
    out = open(outputFasta, "w")
    for (titleStr, seqStr, qualityStr) in FastqIterator(inputFastq):
        out.write(">%s\n%s\n" % (titleStr, seqStr))


def getEdges(start, blockSizes, blockStarts):
    """Returns the left and right edge of a junction.  There must be exactly
    one junction in the line or an exception is raised."""
    if blockSizes.endswith(","):
        sizes = [int(x) for x in blockSizes.split(",")[:-1]]
        starts = [int(x) for x in blockStarts.split(",")[:-1]]
    else:
        sizes = [int(x) for x in blockSizes.split(",")]
        starts = [int(x) for x in blockStarts.split(",")]

    if len(starts) > 2:
        raise Exception("ERROR! one junction per line")

    leftEdge = start + starts[0] + sizes[0]
    rightEdge = start + starts[1]

    return leftEdge, rightEdge  

def divideReads(inputFastq, outputFastq, solexaFormat=False):

    out = open(outputFastq, "w")

    # put the other half of the sequence/quality in the fasta title so we don't have to look it up later
    # (with big fastq files, we may not even be able to put it all in memory to look it up!)
    for (titleStr, seqStr, qualityStr) in FastqIterator(inputFastq):
        # we are going to skip reads with "N" later, may as well not process them through bowtie, etc.
        if seqStr.find("N") >= 0:
            continue
        if solexaFormat:
            qualityStr = convertToQualityStr( convertQualityStr(qualityStr, 64) )
        titleStr = titleStr.replace("|", "_")
        le = len(seqStr)/2
        out.write("@%s|First|%s|%s\n%s\n+\n%s\n" % (titleStr, seqStr[le:], qualityStr[le:], seqStr[:le], qualityStr[:le])) 
        out.write("@%s|Second|%s|%s\n%s\n+\n%s\n" % (titleStr, seqStr[:le], qualityStr[:le], seqStr[le:], qualityStr[le:]))

def fastaDictionary(inFile, chrName=None):
    """return a dictionary of fasta title to sequence strings
    """

    d = {}
    for (title, seq) in FastaIterator(inFile):
        title = title.split()[0]
        if not chrName:
            d[title] = seq
        elif chrName == title:
            d[title] = seq
            return d

    if chrName:
        print "NOT ABLE TO FIND!", chrName
    return d

def trimFastq(infile, outfile, trimLen):
    out = open(outfile, "w")

    for (title, seq, qual) in FastqIterator(infile):
        out.write("@%s\n%s\n+%s\n%s\n" % (title, seq[:trimLen], title, qual[:trimLen]))

    out.close()


def filterFasta(inputFasta, outputFasta, titleStr):

    out = open(outputFasta, "w")

    writing = False
    for line in open(inputFasta):
        if line.startswith(">"):
            if line.find(titleStr) > 0:
                writing = True
                out.write(line)
            else:
                writing = False
        else:
            if writing:
                out.write(line)

#
# Iterator
#
def FastaIterator(fh):
    """return an iterator of Records found in file handle, fh.
    """
    def readTotitle(fh):
        """returns a tuple ([lines before the next title line], next tile line)
        """
        preLines = []
        while True:
            l = fh.readline().strip()
            if l.startswith('>'):
                return (preLines,l)
            elif l == '':
                return preLines,None
            else:
                preLines.append(l)


    if type(fh) in StringTypes:
        fh = file(fh)

    preLines,nextTitleLine =readTotitle(fh)

    while nextTitleLine != None:
        title = nextTitleLine[1:].rstrip()
        preLines,nextTitleLine=readTotitle(fh)
        yield (title,''.join(map(lambda x: x.rstrip(),preLines)))


def FastqIterator(fh):
    """return an iterator of Records found in file handle, fh.
    Original code from Kael and Dale.
    """
    def readTotitle(fh, titleChar):
        """returns a tuple ([lines before the next title line], next tile line)
        """
        preLines = []
        while True:
            l = fh.readline().strip()
            if l.startswith(titleChar):
                return (preLines,l)
            elif l == '':
                return preLines,None
            else:
                preLines.append(l)

    if type(fh) in StringTypes:
        fh = file(fh)

    preLines,nextTitleLine =readTotitle(fh,'@')

    while nextTitleLine != None:
        seqTitle = nextTitleLine[1:].rstrip()
        preLines,nextTitleLine=readTotitle(fh,'+')
        qualTitle = nextTitleLine[1:].rstrip()
        if len(qualTitle.strip()) > 0 and seqTitle != qualTitle:
            print seqTitle
            print preLines
            print qualTitle
            raise hmmErrors.InvalidFastq, "Error in parsing: @title sequence entry must be immediately followed by corresponding +title quality entry."
        seqLines = preLines
        qualLines = []
        for i in range(len(seqLines)): # Quality characters should be the same length as the sequence
            qualLines.append( fh.readline().strip() )

        preLines,nextTitleLine=readTotitle(fh,'@')

        yield (seqTitle, ''.join(seqLines), ''.join(qualLines))


##################################################################################################
#  Utilities to work with output files
##################################################################################################


def filterBedFile(inputBed, outputBed, scoreFilterSingle, scoreFilterMultiple, newName=""):
    """Refilters a set of junction results based on the score.
    inputBed: Should be the full (non-score-filtered) junction results
    outputBed: The filter results
    multipleScore : the score for junctions with support from two or more reads
    singletonScore : the score for junctions with only a single read support
    """

    out = open(outputBed, "w")

    count = 0
    for line in open(inputBed):
        count += 1
        if line.startswith("track"):
            if count > 1:
                continue

            if newName != "":
                pieces = line.split()
                pieces.pop(1)
                pieces.insert(1, "name='%s'" % newName)
                pieces.pop(2)
                pieces.insert(2, "description='%s'" % newName) 
            newTrack = " ".join(pieces)
            out.write(newTrack)
            out.write("\n")
            continue

        pieces = line.split("\t")

        numCollapsed = 0
        if pieces[3].find("|junc=") > 0:
            numCollapsed = int(pieces[3].split("junc=")[-1])

        score = float(pieces[4])
        if (numCollapsed < 2) and (score <= scoreFilterSingle):
            continue
        elif score <= scoreFilterMultiple:
            continue

        out.write("\t".join(pieces))
        # if split on "\t" then "\n" still there.  otherwise need this.
        #out.write("\n")

def getCoveragePerBpFixed(wigFileName):
    """Given a wiggle file name, creates a dictionary with the data in the form
    dictionary[contig][position] = coverage at that position.
    Assumes the wiggle file is in fixedStep format (one column per data line)
    """

    currentContig = ""
    currentPosition = 0
    d = {}

    for line in open(wigFileName):
        if line.startswith("track") or line.startswith("browser"):
            continue

        if line.startswith("fixedStep"):
            currentContig = line.split()[1].split("=")[1]
            currentPosition = 1
            d[currentContig] = []

        else:
            d[currentContig].append(int(line))
            currentPosition = currentPosition + 1

    return d

def getCoveragePerBpVariable(wigFileName):
    """Given a wiggle file name, creates a dictionary with the data in the form
    dictionary[contig][position] = coverage at that position.
    Assumes the wiggle file is in variableStep format (two columns per data line).
    """

    d = {}

    for line in open(wigFileName):
        if line.startswith("track") or line.startswith("browser"):
            continue

        [chr, start, stop, level] = line.split()
        chr = chr.split("|")[1].replace("MAL", "chr")
        level = int(level)
        start = int(start)
        stop = int(stop)

        if not d.has_key(chr):
            d[chr] = {}

        for i in range(start, stop+1):
            d[chr][i] = level

    return d

def getCoverageLevel(inputBed, inputWig, outputTxt, wigTypeIsVariable=True):
    """For each gene (entry) in the bed file, the coverage level in the wig file
    is determined and written to the output."""

    if wigTypeIsVariable:
        cpbp = getCoveragePerBp2(inputWig)
    else:
        cpbp = getCoveragePerBpFixed(inputWig)

    out = open(outputTxt, "w")
    for line in open(inputBed):
        [chr, startStr, stopStr] = line.split()[:3]
        start = int(startStr)
        stop = int(stopStr)

        initialCov = cpbp[chr][start]
        endCov = cpbp[chr][stop]
        covList = []
        for i in range(start, stop):
            covList.append(cpbp[chr][i])
        avgCov = (initialCov + endCov) / 2.0

        out.write("\t".join([chr, startStr, stopStr, str(avgCov), str(max(covList))]))
        out.write("\n")


##################################################################################################
#  Other Utilities
##################################################################################################

def getMinMaxIntron(inputBed, percentBetweenLow, percentBetweenHigh):
    """Analyzes a genome annotation to determine basic intron information. 
    Prints the size of the largest and smallest intron, as well as the total number 
    of introns and the number of introns between the 'percentBetweenLow' and 'percentBetweenHigh'
    paramters.  Also prints the average intron size, average exon size, and the 
    average number of introns per gene.
    """

    minIntron = 1000000
    maxIntron = 0
    countTotal = 0
    countBetween = 0
    sizeList = []
    exonSizeList = []
    numGenes = 0
    for line in open(inputBed):
        if line.startswith("track"):
            continue

        numGenes += 1
        pieces = line.split()
        sizes = [int(x) for x in pieces[10].split(",")[:-1]]
        exonSizeList.extend(sizes)
        starts = [int(x) for x in pieces[11].split(",")[:-1]]
        for i in range(1, len(starts)):
            countTotal += 1
            intronSize = starts[i] - starts[i-1] - sizes[i-1]
            sizeList.append(intronSize)
            if intronSize >= percentBetweenLow and intronSize <= percentBetweenHigh:
                countBetween += 1
            minIntron = min(minIntron, intronSize)
            maxIntron = max(maxIntron, intronSize)

    print "The largest intron was %s and the smallest intron was %s" % (maxIntron, minIntron)
    print "There were %s introns total in %s genes, and %s (%s%%) were between %s and %s (inclusive)" % (countTotal, numGenes, countBetween, 
                                                                                                         (countBetween*100.0/countTotal), percentBetweenLow, percentBetweenHigh)
    print "The average intron size is %.2f" % ( float(sum(sizeList))/len(sizeList))
    print "Average number of introns per gene is %.2f" % ( float(countTotal) / numGenes) 

    print "Average exon size is %.2f" % ( float(sum(exonSizeList)) / len(exonSizeList))


def getFastaFromBed(inputBed, inputGenomeFasta, outputGeneFasta, chrToUse=None):
    """Given a genome annotation and the associated genome fasta file, this function
    creates as fasta file with the sequences of the genes in the input bed file.
    If chrToUse is specified then only the genes on chrToUse are saved."""
    gd = {}
    for (title, seq) in FastaIterator(inputGenomeFasta):
        #rd[title] = seq
        print title
        if chrToUse and title == chrToUse:
            gd[chrToUse] = seq
            #print "     used"
            break

    out = open(outputGeneFasta, "w")
    for line in open(inputBed):
        if line.startswith("track"):
            continue

        pieces = line.split()
        chr = pieces[0]
        start = int(pieces[1])
        sizes = [int(x) for x in pieces[10].split(",")[:-1]]
        starts = [int(x) for x in pieces[11].split(",")[:-1]]
        seq = ""
        for i in range(len(starts)):
            seq = seq + gd[chr][start+starts[i]:start+starts[i]+sizes[i]]
        out.write(">%s\n" % (pieces[3]))
        for i in range(0, len(seq), 60):
            out.write(seq[i:i+60])
            out.write("\n")

def graphQualityPerPosition(inputFastq):
    """Creates a histogram counting the quality value per position in the input fastq file."""

    histD = {}

    count = 0
    for (titleStr, seqStr, qualityStr) in FastqIterator(inputFastq):
        count += 1
        if count < 200000:
            continue
        if count > 1200000:
            break

        qInts = convertQualityStr(qualityStr) 
        for i in range(len(qInts)):
            q = qInts[i]
            if q < 0 or q > 40:
                raise Exception("Invalid quality value %s at position %s of %s" % (q, i, qualityStr))

            if not histD.has_key(i):
                histD[i] = [0]*41

            histD[i][q] += 1

    print "Histogram of quality score per position"
    allk = histD.keys()
    allk.sort()
    for k in allk:
        print "%s|" % k, "|".join(str(x) for x in histD[k])



def bowtieToWig(bowtieFile, trackName, trackDescription, outputName):

    print "Reading file in..."
    allCounts = _readBowtieInput(bowtieFile)
    print allCounts.keys()
    _writeWiggleVar(trackName, trackDescription, allCounts, outputName)
    #_writeWiggle(trackName, trackDescription, allCounts, outputName, win=3)
    #_writeWiggle(trackName, trackDescription, allCounts, outputName)        


def _readBowtieInput(inputFile):
    allCounts = {}
    count = 0
    countIgnored = 0
    for line in open(inputFile):
        count += 1
        if count % 500000 == 0:
            print count

        # save up the number                                                                                    
        [read, strand, chrom, start, seq] = line.split("\t")[:5]         
        start = int(start)
        end = start + len(seq)  

        if not allCounts.has_key(chrom):
            allCounts[chrom] = {}


        for i in range(start, end):
            if allCounts[chrom].has_key(i):
                allCounts[chrom][i] += 1
            else:
                allCounts[chrom][i] = 1

    return allCounts

def _writeWiggle(trackName, trackDescription, allCounts, wigOut, win=1):
    """Writes the allCounts dictionary out to a wiggle file."""
    wigFile = open(wigOut, "w")
    wigFile.write("track type=wiggle_0 name='%s' description='%s' visibility=2\n" % (trackName,
                                                                                     trackDescription))
    for name in allCounts.keys():
        start = 0
        end = max(allCounts[name].keys())

        wigFile.write("fixedStep chrom=%s start=%s step=%s span=%s\n" % (name, start+1, win, win))

        for i in range(start, end):
            if allCounts[name].has_key(i):
                curValue = allCounts[name][i]
            else:
                curValue = 0
            wigFile.write("%s\n" % (curValue))

    wigFile.close()

def _writeWiggleWin(trackName, trackDescription, allCounts, wigOut, win=1):
    """Writes the allCounts dictionary out to a wiggle file."""
    wigFile = open(wigOut, "w")
    wigFile.write("track type=wiggle_0 name='%s' description='%s' visibility=2\n" % (trackName,
                                                                                     trackDescription))

    assert type(win) == type(1), \
           'Error: win must be integer. You gave: %s' % (win)

    for name in allCounts.keys():
        start = 0
        end = max(allCounts[name].keys())

        wigFile.write("fixedStep chrom=%s start=%s step=%s span=%s\n" % (name, start+1, win, win))
        i=0

        while i <= end:
            winCounts = []
            for pos in range(i,i+win):
                if allCounts[name].has_key(pos):
                    winCounts.append(allCounts[name][pos])
                else:
                    winCounts.append(0)
            winAvg = sum(winCounts)/float(len(winCounts))
            wigFile.write("%.4f\n" % (winAvg))
            i+=win

    wigFile.close()

def _writeWiggleVar(trackName, trackDescription, allCounts, wigOut):
    """Writes the allCounts dictionary out to a wiggle file."""
    wigFile = open(wigOut, "w")
    wigFile.write("track type=wiggle_0 name='%s' description='%s' visibility=2\n" % (trackName,
                                                                                     trackDescription))

    for name in allCounts.keys():
        start = 0
        end = max(allCounts[name].keys())

        wigFile.write("variableStep chrom=%s span=1\n" % (name))

        for pos in sorted(allCounts[name].keys()):
            wigFile.write("%s\t%s\n" % (pos,allCounts[name][pos]))


    wigFile.close()    


def readJunctionsFromBed(inputBed, saveWholeLine=False, wiggle=0):
    junct = {}
    for line in open(inputBed):
        if line.startswith("track"):
            continue

        if len(line) < 2:
            continue

        [chr, start, stop, name, score, strand, thStart, thStop, rgb, blockCount, blockSizes, blockStarts] = line.split()
        if not junct.has_key(chr):
            junct[chr] = {}

        start = int(start)
        stop = int(stop)
        blockCount = int(blockCount)
        if blockSizes.endswith(","):
            sizes = [int(x) for x in blockSizes.split(",")[:-1]]
            starts = [int(x) for x in blockStarts.split(",")[:-1]]
        else:
            sizes = [int(x) for x in blockSizes.split(",")]
            starts = [int(x) for x in blockStarts.split(",")]
        for i in range(1, blockCount):
            leftEdge = start + starts[i-1] + sizes[i-1]
            rightEdge = start + starts[i]  
            if saveWholeLine:
                if not junct[chr].has_key((leftEdge, rightEdge)):
                    junct[chr][(leftEdge, rightEdge)] = []
                junct[chr][(leftEdge, rightEdge)].append(line.strip())
            else:
                if rightEdge - leftEdge > 15:
                    found = False
                    for wL in range(-wiggle,wiggle):
                        for wR in range(-wiggle,wiggle):
                            if junct[chr].has_key( (leftEdge+wL, rightEdge+wR)):
                                found = True
                                break
                    if not found:
                        junct[chr][ (leftEdge, rightEdge) ] = (name, score, strand)

    return junct

def hasJunction(junc, chr, leftEdge, rightEdge, wiggle=0):
    """Determines whether the given junction dictionary has a junction.
    wiggle is the amount plus/minus to count as the 'same place'."""

    for i in range(leftEdge-wiggle, leftEdge+wiggle+1):
        for j in range(rightEdge-wiggle, rightEdge+wiggle+1):
            try:
                if junc[chr].has_key( (i, j) ):
                    return True
            except KeyError:
                return False

    return False


def measureSpecificity(foundJunctionsBed, estJunctionsBed, wiggle=0, goodFile=None, badFile=None):

    # read in all the est junctions
    ests = readJunctionsFromBed(estJunctionsBed)

    return measureSpecificityUsingDict(foundJunctionsBed, ests, wiggle, goodFile, badFile)

def measureSpecificityUsingDict(foundJunctionsBed, ests, wiggle=0, goodFile=None, badFile=None):
    """Prints the percent of found junctions that overlap with ests."""

    overlaps = 0
    noOverlaps = 0

    # do it slow for now
    # read in all the junctions we found
    found = {}  # found[chr][(left, right)] = count
    for line in open(foundJunctionsBed):
        if line.startswith("track"):
            continue

        if len(line) < 3:
            continue

        pieces = line.split("\t")
        #print pieces
        if len(pieces) == 0:
            continue

        if pieces[0].startswith("Plasmodium_falciparum"):
            pieces[0] = pieces[0].split("|")[1].replace("MAL", "chr")

        if pieces[0].startswith("psu|Pf"):
            pieces[0] = "chr" + str(int(pieces[0].split()[0].split("|")[1].split("_")[1]))


        if not found.has_key(pieces[0]):
            found[pieces[0]] = {}

        leftEdge, rightEdge = getEdges(int(pieces[1]), pieces[10], pieces[11])
        if not found[pieces[0]].has_key( (leftEdge,rightEdge) ):
            found[pieces[0]][(leftEdge,rightEdge)] = 0
        found[pieces[0]][(leftEdge,rightEdge)] += 1

    # for every one of our junction, do they overlap with an est?
    if goodFile != None and badFile != None:
        goodOut = open(goodFile, "w")
        badOut = open(badFile, "w")

    for chr, edgeDict in found.iteritems():
        for (leftEdge, rightEdge) in edgeDict.keys():
            foundOne = False
            for x in range(leftEdge-wiggle, leftEdge+wiggle+1):
                for y in range(rightEdge-wiggle, rightEdge+wiggle+1):
                    if ests.has_key(chr):
                        if ests[chr].has_key( (x, y) ):
                            foundOne = True
            if foundOne:
                overlaps += found[chr][(leftEdge,rightEdge)]
                if goodFile != None:
                    goodOut.write("\t".join([chr, str(leftEdge), str(rightEdge), str(overlaps)]))
                    goodOut.write("\n")
            else:
                noOverlaps += found[chr][(leftEdge,rightEdge)]
                if badFile != None:
                    badOut.write("\t".join([chr, str(leftEdge), str(rightEdge), str(noOverlaps)]))
                    badOut.write("\n")

    if (noOverlaps + overlaps) > 0:
        print "%s overlapped but %s did not.  %.2d%% overlapped" % (overlaps, noOverlaps, (overlaps*100.0)/(noOverlaps+overlaps))

        return (overlaps*100.0) / (noOverlaps + overlaps)
    else:
        print "No junctions found!"
        return 0

def generateROCcurve(junctionBed, ests, stepSize=100, wiggle=0):
    """Generates an ROC curve to measure good/bad junctions in the junction bed file.

    The ROC curve is calculated from the 'good' and 'bad' matches, with the criteria being the score cut-off.
    """

    goodNames, badNames = findGoodBadReads(junctionBed, ests, wiggle)

    rocTable = []
    scores = []
    scoredDict = {}
    countGood = 0
    countBad = 0
    countStacked = 0
    countUnclear = 0

    for line in open(junctionBed):
        if line.startswith("track"):
            continue
        [chr, start, stop, name, score] = line.split()[:5]
        score = float(score)
        if goodNames.has_key(name):
            isGood = 1 #True
        elif badNames.has_key(name):
            isGood = 0 #False
        else:
            continue

        scores.append(score)

        if isGood:
            countGood += 1

        else:
            countBad += 1

        scoredDict[name] = (score, isGood)

    print "Counted %s good and %s bad" % (countGood, countBad)

    #print scores[:3]
    scores.sort(reverse=True)
    #print scores[:5]
    #use every 100th score
    for x in scores[::stepSize]:
        numGoodAbove = 0
        numBadAbove = 0
        for k, (score, real) in scoredDict.iteritems():
            if score >= x:
                #print real, score, k, pairsDict[k]
                if real:
                    numGoodAbove += 1
                else:
                    numBadAbove += 1

        #print x, numGoodAbove, numBadAbove, numGoodAbove / float(countGood) , numBadAbove / float(countBad)
        rocTable.append((x, numGoodAbove, numBadAbove, numGoodAbove / float(countGood) , numBadAbove / float(countBad)))
    return rocTable,scoredDict

def findGoodBadReads(junctionBed, ests, wiggle):
    """Finds 'good' and 'bad' reads, only using the estBed as a guide."""

    ests = readJunctionsFromBed(ests)

    goodNames = {}
    badNames = {}
    for line in open(junctionBed):
        if line.startswith("track"):
            continue

        [chr, start, stop, name, score, strand, thStart, thStop, rgb, blockCount, blockSizes, blockStarts] = line.split()

        hasEst = testIfEst(ests, chr, int(start), int(blockCount), blockSizes, blockStarts, wiggle)

        if hasEst:
            goodNames[name] = 1
        else:
            badNames[name] = 1

    return goodNames, badNames

def vennDiagram(bed1File, bed2File, only1Output=None, only2Output=None, bothOutput=None):
    """Counts the number of splice junctions unique to bed1, unique to bed2, and shared in both."""

    bed1 = readJunctionsFromBed(bed1File, True)
    bed2 = readJunctionsFromBed(bed2File, True)

    count1 = 0
    count2 = 0
    countBoth = 0

    out1 = None
    if only1Output:
        out1 = open(only1Output, "w")
    out2 = None
    if only2Output:
        out2 = open(only2Output, "w")
    both = None
    if bothOutput:
        both = open(bothOutput, "w")

    for chr, chrJunct in bed1.iteritems():
        for (start,stop) in chrJunct:
            if bed2.has_key(chr):
                if bed2[chr].has_key( (start, stop) ):
                    if both:
                        for line in bed1[chr][(start,stop)]:
                            both.write(line)
                            both.write("\n")
                    del bed2[chr][(start,stop)]
                    countBoth += 1
                else:
                    count1 += 1
                    if out1:
                        line = bed1[chr][(start,stop)][0]
                        pieces = line.split()
                        bedVals = [chr, start-10, stop+10, pieces[3], pieces[4], pieces[5], start-10, stop+10, pieces[8], pieces[9],
                                   "10,10", "0,%s"%(stop-start+10)]
                        out1.write("\t".join(str(x) for x in bedVals))
                        out1.write("\n")
                        #for line in bed1[chr][(start, stop)]:
                        #    out1.write(line)
                        #    out1.write("\n")
            else:
                count1 += 1
                if out1:
                    line = bed1[chr][(start,stop)][0]
                    pieces = line.split()
                    bedVals = [chr, start-10, stop+10, pieces[3], pieces[4], pieces[5], start-10, stop+10, pieces[8], "2",
                               "10,10", "0,%s"%(stop-start+10)]
                    out1.write("\t".join(str(x) for x in bedVals))
                    out1.write("\n")
                    #for line in bed1[chr][(start, stop)]:
                    #    out1.write(line)
                    #    out1.write("\n")

    #print
    #print
    #print

    count2 = sum( len(chrJunct) for chrJunct in bed2.values())
    if out2:
        for chr, chrJunct in bed2.iteritems():
            for (start,stop) in chrJunct:
                line = bed2[chr][(start,stop)][0]
                pieces = line.split()
                bedVals = [chr, start-10, stop+10, pieces[3], pieces[4], pieces[5], start-10, stop+10, pieces[8], "2",
                           "10,10", "0,%s"%(stop-start+10)]
                out2.write("\t".join(str(x) for x in bedVals))
                out2.write("\n")
                #for line in bed2[chr][(start, stop)]:
                #    out2.write(line)
                #    out2.write("\n")

    print "There were %s in both, %s in the first one and %s in the second one" % (countBoth, count1, count2)



def convertAedesContigs(genomeFasta, junctionBed, fixedName):

    # read all the CH names and their corresponding supercontigs into a dictionary
    d = {}
    for (title, seq) in FastaIterator(genomeFasta):
        oldName = title.split()[0]
        newName = title.split()[-1]
        d[oldName] = newName   

    out = open(fixedName, "w")
    # go through junction bed and replace
    for line in open(junctionBed):
        if line.startswith("track"):
            out.write(line)
            continue

        pieces = line.split()
        #print "Looking for %s" % (pieces[0])
        pieces[0] = d[pieces[0]]
        out.write("\t".join(pieces))
        out.write("\n")

#def getEdges(start, blockSizes, blockStarts):
    #"""Returns the left and right edge of a junction.  There must be exactly
    #one junction in the line or an exception is raised."""
    #if blockSizes.endswith(","):
        #sizes = [int(x) for x in blockSizes.split(",")[:-1]]
        #starts = [int(x) for x in blockStarts.split(",")[:-1]]
    #else:
        #sizes = [int(x) for x in blockSizes.split(",")]
        #starts = [int(x) for x in blockStarts.split(",")]

    #if len(starts) > 2:
        #raise Exception("ERROR! one junction per line")

    #leftEdge = start + starts[0] + sizes[0]
    #rightEdge = start + starts[1]

    #return leftEdge, rightEdge 

def getEdges(start, blockSizes, blockStarts):
    """Returns the left and right edge of a junction.  There must be exactly
    one junction in the line or an exception is raised.

    The start, blockSizes, and blockStarts are from BED file format."""
    if blockSizes.endswith(","):
        sizes = [int(x) for x in blockSizes.split(",")[:-1]]
        starts = [int(x) for x in blockStarts.split(",")[:-1]]
    else:
        sizes = [int(x) for x in blockSizes.split(",")]
        starts = [int(x) for x in blockStarts.split(",")]

    if len(starts) > 2:
        raise hmmErrors.InvalidInputException("getEdges: ERROR! onejunction per line.  Input was: %s, %s, %s" % (start, blockSizes,
blockStarts))

    leftEdge = int(start) + starts[0] + sizes[0]
    rightEdge = int(start) + starts[1]

    return leftEdge, rightEdge



############ Added via "May 26, 2010 8:49:08 AM PDT" email ############

def testIfEst(ests, chr, start, blockCount, blockSizes, blockStarts, wiggle=0):
    """Determines whether the ests dictionary contains the junction described by the                                                        
    chr, start, blockCount, blockSizes."""
    if int(blockCount) != 2:
        print "ERROR!  the block count isn't 2!  %s, %s, %s, %s, %s" % (chr, start, blockCount, blockSizes, blockStarts)
        return False

    (leftEdge, rightEdge) = getEdges(start, blockSizes, blockStarts)

    return hasJunction(ests, chr, leftEdge, rightEdge, wiggle)


def collapseCloseJunctions(inputBed, outputBed, withinBp):
    """Collapses junctions that are within 'withinBp' bp of each other.  The most probable junction                                         
    of the close junctions are used.  Probabilities are converted appropriately (see collapseJunctions                                      
    for description).                                                                                                                       
    Also collapses identical junctions in the same step.                                                                                    
    Use withinBp=0 to only collapse identical junctions                                                                                     
    """

    junct = _readAndCombine(inputBed, withinBp)

    # then go through each dictionary item                                                                                                  
    #   -- if the list is length 1, just write the line                                                                                     
    #   -- if the list is longer than 1, combine probabilities, adjust the name to include the number, and write 1                          
    out = open(outputBed, "w")
    out.write("track name=collapsedJunctions description='Junctions' useScore=1\n")
    for chr, junctions in junct.iteritems():
        for (leftEdge, x,rightEdge, y, intronLength), junctionList in junctions.iteritems():
            if len(junctionList) == 1:
                out.write(junctionList[0][1].strip())
                #pieces = junctionList[0][0].split()                                                                                        
                #pieces.pop(4)                                                                                                              
                #pieces.insert(4, "100")                                                                                                    
                #out.write("\t".join(pieces))                                                                                               
                out.write("\n")
            else:
                out.write(_combineLines(junctionList, leftEdge, rightEdge))
                out.write("\n")


def _readAndCombine(inputBed, withinBp):
    """Helper for collapseCloseJunctions reads in the inputBed into a dictionary."""
    junct = {}

    # collapse a                                                                                                                            
    count = 0
    for line in open(inputBed):
        count += 1
        #if count % 100000==0:                                                                                                              
        #    print count                                                                                                                    
        if line.startswith("track"):
            #out.write(line.strip())                                                                                                        
            #out.write(" useScore=1\n")                                                                                                     
            continue

        [chr, start, stop, name, score, strand, thStart, thStop, rgb, blockCount, blockSizes, blockStarts] = line.split("\t")
        score = float(score)
        if not junct.has_key(chr):
            junct[chr] = {}

        if int(blockCount) != 2:
            #print "Illegal line does not have 2 blocks"                                                                                    
            #print line                                                                                                                     
            continue

        start = int(start)
        stop = int(stop)
        [size1, size2] = [int(x) for x in blockSizes.split(",")[:2]]
        [start1, start2] = [int(x) for x in blockStarts.split(",")[:2]]
        leftEdge = start + size1
        rightEdge = start + start2  # start2 is relative to chr start                                                                       
        intronLength = rightEdge - leftEdge

        toCombine = []
        for (other) in junct[chr].keys():
            (otherMinLeft, otherMaxLeft, otherMinRight, otherMaxRight, otherLength) = other
            if otherLength != intronLength:
                continue

            if otherMaxLeft < (leftEdge-withinBp) or otherMinLeft > (leftEdge+withinBp):
                continue

            if otherMaxRight < (rightEdge-withinBp) or otherMinRight > (rightEdge+withinBp):
                continue

            toCombine.append(other)

        allLines = [ (score, line, leftEdge, rightEdge) ]
        minLeft = maxLeft = leftEdge
        minRight = maxRight = rightEdge
        for (other) in toCombine:
            (otherMinLeft, otherMaxLeft, otherMinRight, otherMaxRight, intronLength) = other
            minLeft = min(minLeft, otherMinLeft)
            maxLeft = max(maxLeft, otherMaxLeft)
            minRight = min(minRight, otherMinRight)
            maxRight = max(maxRight, otherMaxRight)

            allLines.extend(junct[chr][other])
            del junct[chr][other]

        junct[chr][ (minLeft, maxLeft, minRight, maxRight, intronLength) ] = allLines

    return junct

def _combineLines(junctionList, leftEdge, rightEdge):
    """Helper method for collapseCloseJunctions.                                                                                            
    Takes several junction lines and combines them into a single line, then returns that line.                                              
    The score for the combined junction is increased accordingly.                                                                           
    """

    # the only things we have to adjust are the score and the name                                                                          
    scoreSoFar = 0
    previousNumBasesCovered = 0
    minStart = -1
    maxStop = -1
    # keep a count of all the left/right edges so we can pick the most common                                                               
    leftEdgeCount = {}
    rightEdgeCount = {}

    # first we need to sort the junctionList by score, in reverse (highest score first)                                                     
    junctionList.sort(reverse=True)
    name = junctionList[0][1].split("\t")[3]
    countPlus = 0
    countMinus = 0
    numJunctions = 0
    for (score, line, leftEdge, rightEdge) in junctionList:
        pieces = line.split("\t")
        #print "previous start: %s, this start: %s" % (minStart, pieces[1])                                                                 
        start = int(pieces[1])
        stop = int(pieces[2])
        strand = pieces[5]
        if strand == "+":
            countPlus += 1
        elif strand == "-":
            countMinus += 1
        else:
            print line
            raise hmmErrors.InvalidInputException("ERROR!  strand value %s not valid. " % strand)

        name = pieces[3]
        if pieces[3].find("|junc=") > 0:
            numCollapsed = int(pieces[3].split("junc=")[-1])
            numJunctions += numCollapsed
        else:
            numJunctions += 1

        if previousNumBasesCovered == 0:
            previousNumBasesCovered = (leftEdge - start) + (stop - rightEdge)
            scoreSoFar = float(pieces[4])
        else:
            # we only want to count bases grown on the outer sides (start and stop) and ignore bases on the inner edges                     
            newBases = max(0, minStart-start) + max(0, stop-maxStop)
            scoreSoFar = scoreSoFar + (newBases / float(newBases+previousNumBasesCovered) ) * float(pieces[4])
            previousNumBasesCovered = previousNumBasesCovered + newBases


        if minStart > 0:
            minStart = min(minStart, start)
        else:
            minStart = start

        maxStop = max(maxStop, stop)

        if not leftEdgeCount.has_key(leftEdge):
            leftEdgeCount[leftEdge] = 0
        leftEdgeCount[leftEdge] += 1
        if not rightEdgeCount.has_key(rightEdge):
            rightEdgeCount[rightEdge] = 0
        rightEdgeCount[rightEdge] += 1

    maxLeft = max(leftEdgeCount.values())
    for k, v in leftEdgeCount.iteritems():
        if v == maxLeft:
            useLeft = k
            break
    maxRight = max(rightEdgeCount.values())
    for k, v in rightEdgeCount.iteritems():
        if v == maxRight:
            useRight = k
            break
    if countPlus >= countMinus:
        strand = "+"
    else:
        strand = "-"

    pieces = junctionList[0][1].split("\t")
    namePieces = pieces[3].split("|")
    rootName = ""
    for piece in namePieces:
        if not piece.startswith("junc="):
            rootName += piece + "|"
    finalName = rootName + ("junc=%s" % numJunctions)
    blockStarts = "0,%s," % (useRight - minStart)
    blockSizes = "%s,%s," % ( (useLeft-minStart), (maxStop-useRight) )

    return "\t".join(str(x) for x in [pieces[0], minStart, maxStop, finalName, scoreSoFar,
                                      strand, minStart, maxStop, "0", "2", blockSizes, blockStarts
                                      ])



########### Below be DRAGONS!! ###############
########### (gus's additions)  ###############












