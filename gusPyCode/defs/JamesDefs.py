import re
import random
import copy
import os
import collections
from gusPyCode.defs import xpermutations



#=========================================================================
# 02/20/10
def detect_1D_overlap(coords1,coords2):
    """Returns TRUE if coords1 overlaps with coords2.
    """
    coords1[0],coords1[1] = int(coords1[0]),int(coords1[1])
    coords2[0],coords2[1] = int(coords2[0]),int(coords2[1])
    # +++ Validate Usr Input +++
    assert (len(coords1)==2) and (len(coords2)==2), \
           "** ERROR: coords1 and coords2 must be lists of length 2! **"
    # +++ Sort Coords +++
    coords1.sort()
    coords2.sort()
    # +++ Classify Coords +++
    if (coords1[1]-coords1[0]) <= (coords2[1]-coords2[0]):
        shorter = coords1
        longer  = coords2
    else:
        shorter = coords2
        longer  = coords1
    # +++  +++
    lEdge = (shorter[0]-longer[0] >= 0) and (shorter[0]-longer[1] <= 0)
    rEdge = (shorter[1]-longer[0] >= 0) and (shorter[1]-longer[1] <= 0)
    # -- did we get a hit? --
    return (lEdge or rEdge)

#=========================================================================
# 01/06/10
def tbl2nmdTup(filePath, name):
    data = map(lambda l: l.strip('\n').split('\t'), open(filePath).readlines())

    if data[0][0].startswith('#'):
        data[0][0] = data[0][0][1:]
        header = data.pop(0)
    else:
        print 'ERROR:  The first line of table file must be the header and begin with a "#".'
        exit(1)
    # ---- Replace ' ' with '_' in header ----
    for i in range(len(header)):
        header[i] = header[i].replace(' ','_')
        
    namedTab = collections.namedtuple(name,' '.join(header))
    return [namedTab._make(x) for x in data]

#=========================================================================
# 12/21/09

def mkdirp(path):
    """Create new dir while creating any parent dirs in the path as needed.
    """
    if not os.path.isdir(path):
        os.makedirs(path)
#=========================================================================
# 11/23/09
# From Titus Brown's gff parser:
class Bag(dict):
    """dict-like class that supports attribute access as well as getitem.

    >>> x = Bag()
    >>> x['foo'] = 'bar'
    >>> x.foo
    'bar'
    
    """
    def __init__(self, *args, **kw):
        dict.__init__(self, *args, **kw)
        for k in self.keys():
            self.__dict__[k] = self.__getitem__(k)

    def __setitem__(self, k, v):
        dict.__setitem__(self, k, v)
        self.__dict__[k] = v

#=========================================================================
# 11/07/09
def readDelimFile(filePath,delimiter='\t'):
    return map(lambda l: l.strip('\n').split(delimiter), open(filePath, 'rU').readlines())

#=========================================================================
# 10/29/09
class DotDict(dict):
    """
    Class to package multiple objects with before pickling to make them callable
    with the "unPickledObj.subObj" convention instead of using a dictionary.
    """
    def __getattr__(self, attr):
        return self.get(attr, None)
    __setattr__= dict.__setitem__
    __delattr__= dict.__delitem__


#=========================================================================
# 10/08/09
def initList(numOfIndices,initVal):
    """
    Returns a list of length len(numOfIndices) initialized with fresh invocations of the
    empty(or initialized) data type supplied with initVal.
    """
    
    rList = []
    for i in range(numOfIndices):
        if (type(initVal) == type([])) or (type(initVal) == type({})) or (type(initVal) == type(set())):
            rList.append(eval(repr(initVal)))  # prevent init'ing all indices with ref to same obj.
        else:
            rList.append(initVal)
    return rList
#=========================================================================

#=========================================================================
# 07/19/09
def odd_or_even(integer):
    assert type(integer) == type(1), 'Error: odd_or_even only takes integers. You gave: %s' % (integer)
    if integer % 2 == 0:
        return "even"
    else:
        return "odd"
#=========================================================================


#=========================================================================
# 06/05/09
#def cmpListOfLists(a,b):
    
#=========================================================================
# 08/26/09
def randFromList_noRplcMulti(data, numToPull):
    assert numToPull <= len(data), \
           'ERROR in randFromList_noRplcMulti: numToPul > len(data)'
    data = data[:]
    
    theChoosenOnes = []
    for i in range(numToPull):
        
        if data != []:
            index = random.randint(0, len(data) - 1)
            elem = data[index]
            data[index] = data[-1]
            del data[-1]
            theChoosenOnes.append(elem)
        else:
            raise
        
    return theChoosenOnes
#=========================================================================
# 06/04/09
def randFromList_noReplace(data):
    if data != []:
        index = random.randint(0, len(data) - 1)
        elem = data[index]
        data[index] = data[-1]
        del data[-1]
        return elem
    else:
        return data
#=========================================================================
# 05/28/09
def combineOrthoTabs(orthoTabsList, infer=False):
    """
    Take list of ortholog def lists in format:list[GeneIDspecies1,GeneIDspecies2]
    Return a list of orthologs found in all species.
    
    If infer: allow inferences of N-way 1:1 orthos even if a members 1:1 pair is missing
    in a genomePair's list.  Else: require that the 1:1 relationship be explicitly defined
    in all genomes.
    """
    
    orthoTabs = copy.deepcopy(orthoTabsList)
    
    # If not infer: convert orthoTabsList to list of sets of sets
    # (since infer defualts to False this is default behavior)
    if infer == False:
        for i in range(len(orthoTabsList)):
            for j in range(len(orthoTabsList[i])):
                orthoTabsList[i][j] = frozenset(orthoTabsList[i][j])
            orthoTabsList[i] = frozenset(orthoTabsList[i]) 
    
    
    # Learn how many Genomes we are dealing with
    # and log the species prefixes
    
    # Get representative geneNames
    geneNames = []
    for i in orthoTabs:
        geneNames.extend(i[0])

    # Extract letter portion of gene tokens
    genomes = set([])
    prefixRegEx = re.compile('^\D+', re.IGNORECASE)
    for name in geneNames:
        assert prefixRegEx.match(name) != None , 'geneNames(%s) do not seem to have form "PREFIX000000"' % (name)
        genomes.add(prefixRegEx.match(name).group())
    genomes = list(genomes)
        
    # Move longest list to first position and use _IT_ as root list for merging
    orthoTabSizes = []
    for i in orthoTabs:
        orthoTabSizes.append(len(i))
    ##print '%s' % (str(orthoTabSizes))
    for i in range(len(orthoTabs)):
        if len(orthoTabs[i]) == max(orthoTabSizes):
            biggest = orthoTabs.pop(i)
            orthoTabs.insert(0,biggest)
    
    orthoTabSizes2 = []
    for i in orthoTabs:
        orthoTabSizes2.append(len(i))       
    ##print '%s' % (str(orthoTabSizes2))
    
    # Set-ify all orthoPairs
    for i in range(len(orthoTabs)):
        for j in range(len(orthoTabs[i])):
            orthoTabs[i][j] = set(orthoTabs[i][j])
    
    # Merge sets if any geneIDs are in common
    for i in range(1,len(orthoTabs)):                          # Start w/ 2nd orthoList and compare it and all after to 1st list
        for j in range(len(orthoTabs[i])):                     # For every j in 2nd or higher 
            for k in range(len(orthoTabs[0])):                 # Test every k in 1st orthoList for j
                if orthoTabs[0][k]&orthoTabs[i][j]:            # If there is intersection bt k and ij,
                    orthoTabs[0][k].update(orthoTabs[i][j])    # merge k and ij
                    
    # Remove redundancy and rename merged list for ease of use
    print 'len of mergedIDs list = %s' % (len(orthoTabs[0]))
    mergedIDs = orthoTabs[0]
    for i in range(len(mergedIDs)):
        mergedIDs[i] = frozenset(mergedIDs[i])
    mergedIDs = list(set(mergedIDs))
    print 'len of nr mergedIDs list = %s' % (len(mergedIDs))
    

    
    
    # Sort each orthoSet
    for i in range(len(mergedIDs)):
        mergedIDs[i] = list(mergedIDs[i])
        mergedIDs[i].sort()
        
    # Retain only those relationships that are the correct length
    # AND do not repeat Species prefixs ('AGAP')    
    filteredIDs = []
    for orthoSet in mergedIDs:
        prefixes = []
        for ID in orthoSet:
            for g in genomes:
                if ID.startswith(g):
                    prefixes.append(g)
                    continue
        if len(prefixes) != len(genomes):  # can get same length without having same constiuents 
            continue
        if len(prefixes) == len(set(prefixes)):
            if infer != False:
                filteredIDs.append(orthoSet)
            else:
                orthoPairCount = 0
                for pair in xpermutations.xuniqueCombinations(orthoSet,2):
                    for genomePairs in orthoTabsList:
                        if set(pair) not in genomePairs: continue
                        else: orthoPairCount+=1
                if orthoPairCount == len(orthoTabsList):
                    filteredIDs.append(orthoSet)
                    
    
        
    # test for repeats in each genome's list of final names
    genomeNameLists = []
    for i in range(len(filteredIDs[0])):
        tempL = []
        for j in range(len(filteredIDs)):
            tempL.append(filteredIDs[j][i])
        genomeNameLists.append(tempL)
    oneGenomeNRs = []    
    for each in genomeNameLists:
        oneGenomeNRs.append(len(set(each)))
    c1 = 0
    for each in oneGenomeNRs:
        c1 += 1
        print '%s NRs in genome %s' % (each,c1)
            
            
    return filteredIDs

#=========================================================================
# 12/13/08
def loadXXmiRNAs(filePath):
    f = open(filePath, 'rU').readlines()
    f = ''.join(f)
    
    miRNAs = f.split('>')
    miRNAs.pop(0)
    
    #split miRNA indexes into their own lists
    for i in range(0,len(miRNAs)):
        miRNAs[i] = miRNAs[i].strip('\n').split('\n')
    
    return miRNAs
        


#=========================================================================
# 12/10/08
def removeCommentLines(listOfLines,commentChar):
    cleansed = []
    for l in listOfLines:
        if l.startswith(commentChar):
            continue
        else:
            cleansed.append(l)
    return cleansed

#=========================================================================
#def overlapRegEx(reObj,seq,skip=1,startPos=0,count=0,posList=[],hitList=[]):
#	""" This def does not work.  It does not seem to kill its vars and the info 
#	gets used later when it should not still exist.  Dont use in current form."""
#    ## This is bc these two lists kept surviving multiple calls to this func
#    ## >There must be better way but I am bored with this.
#    #posList = eval(str(posList))
#    #hitList = eval(str(hitList))
#    h = reObj.search(seq,startPos)
#    print 'matchObj:%s' %  (h)
#    ##if h:
#        ##print '%s\t%s\t%s' % (h.start(),h.end(),h.group())
#        ##x=1
#    if h:
#        posList.append(h.start())
#        hitList.append(h.group())
#        count += 1
#        startPos += h.start()+skip
#        #posList = str(posList)
#        #hitList = str(hitList)
#        returnList = overlapRegEx(reObj,seq,skip,startPos,count,posList,hitList)
#        return returnList
#    posList = str(posList)
#    hitList = str(hitList)
#    return str([count,eval(posList),eval(hitList)])
#
#=========================================================================
def cpu():
    import resource
    return (resource.getrusage(resource.RUSAGE_SELF).ru_utime+
            resource.getrusage(resource.RUSAGE_SELF).ru_stime)

#=========================================================================



def nrListBySets(importList):
    """
    takes list
    returns [int(num of items removed), nrList]
    """
    from sets import Set
    startLen = len(importList)
    nrList = list(Set(importList))
    endLen = len(nrList)
    return [startLen-endLen, nrList]




def fastaFileToBioSeqDict(pathToFastaFile, Alphabet='IUPACAmbiguousDNA', splitOn=None, indexAsName=0):
    from Bio.Alphabet.IUPAC import Alphabet
    from Bio import SeqIO
    #  Populate a Dict with Seq objs for Anopheles boundary seqs
    #  What follows directly is a klugde to get all my seqDict vals to have an IUPAC alphabet
    listOfSeqs = list(SeqIO.parse(open(pathToFastaFile, "rU"), "fasta"))
    for record in listOfSeqs :
        record.seq.alphabet = Alphabet

    dictOfSeqs = SeqIO.to_dict(listOfSeqs, key_function = lambda rec : rec.description.split(splitOn)[indexAsName])

    return dictOfSeqs

#=========================================================================
compl_iupacdict = {'A':'T',
                   'C':'G',
                   'G':'C',
                   'T':'A',
                   'M':'K',
                   'R':'Y',
                   'W':'W',
                   'S':'S',
                   'Y':'R',
                   'K':'M',
                   'V':'B',
                   'H':'D',
                   'D':'H',
                   'B':'V',
                   'X':'X',
                   'N':'N',
                   '-':'-'}

def compliment(motif, compl_iupacdict):
    compl_motif = ""
    for i in range(0,len(motif)):
        letter = motif[i]
        compl_motif = compl_motif + compl_iupacdict[letter]
    return compl_motif

def reverse(text):
    return text[::-1]

def revComp(seq):
    revCompSeq = reverse(compliment(seq, compl_iupacdict))
    return revCompSeq
#=========================================================================

def iupacList_2_regExList(motifList):
    i = 0
    while i < len(motifList):
        motifList[i] = [motifList[i], iupac2regex(motifList[i])]
        i += 1

def iupac2fwdRevCmpTRegExObj(motif):
    """Takes motif in IUPAC form (WGATAR).  Returns a compiled python regular expression
    object recognizing either upper or lowercase of the fwd or revComp version of the motif."""
    
    import re
    
    # list of few and rev
    fwdAndRevCp = [motif,revComp(motif)]
    
    # convert WGATAR to [AT]GATA[AG] for fwd and rev
    for i in range(2):
        fwdAndRevCp[i] = iupac2regex(fwdAndRevCp[i])
    
    # create string to feed to re
    reString = '(%s|%s)' % (fwdAndRevCp[0],fwdAndRevCp[1])
    
    # return compiled regEx Obj
    return re.compile(reString,re.IGNORECASE)


def iupac2regex(motif):
    """Convert 'WGATAR' to '[AT]GATA[AG]'"""
    
    iupacdict = {'A':'A',
                 'C':'C',
                 'G':'G',
                 'T':'T',
                 'M':'[AC]',
                 'R':'[AG]',
                 'W':'[AT]',
                 'S':'[CG]',
                 'Y':'[CT]',
                 'K':'[GT]',
                 'V':'[ACG]',
                 'H':'[ACT]',
                 'D':'[AGT]',
                 'B':'[CGT]',
                 'X':'[ACGT]',
                 'N':'[ACGT]'}

    transl_motif = ""
    for i in range(0,len(motif)):
        letter = motif[i]
        transl_motif = transl_motif + iupacdict[letter]
    return transl_motif

#=========================================================================

def reOrderDelimitedList(listOfDelimStrings, Delimiter, newOrderList):
    """
    Takes:	- one list delimited with a unique char and a list with new field indexes:
	                reOrderDelimitedList(listOfDelimStrings, strToSplitOn, newOrderList)
			- NOTE: any trailing newlines should be REMOVED

    Does:	- explodes delimited list item into a list of its fields IN PLACE with explodeDelimitedList()
			- uses a while loop to build new STRING with desired field order using original delimiter
			- appends this string to returnList

    Returns:		- new list of delimited strings with reorganized fields
    """

    returnList = []
    #  Convert list into list of lists splitting on delimiter
    explodeDelimitedList(listOfDelimStrings, Delimiter)

    for each in listOfDelimStrings:

        #  Build new string based on order of indexes given in newOrderList
        newDelimitedStr = ''
        i = 0
        while i < len(newOrderList):
            newDelimitedStr = newDelimitedStr+each[newOrderList[i]]+Delimiter
            i = i+1
        #  Remove trailing delimiter and append to returnList
        newDelimitedStr.rstrip(Delimiter)
        returnList.append(newDelimitedStr)

    return returnList

#=========================================================================

def explodeDelimitedList(listOfDelimStrings, Delimiter):
    """
    Takes:		- one list delimited with a unique char: explodeListOfDelimStr(list, 'char')
                        - NOTE: any trailing newlines should be REMOVED

    Does:		- explodes delimited list item into a list of its fields IN PLACE

    Returns:		- Nothing.  It acts on the source list in place.
    """

    listLen = len(listOfDelimStrings)
    i = 0
    while i < listLen:
        listOfDelimStrings[i] = listOfDelimStrings[i].rstrip('\n')
        listOfDelimStrings[i] = listOfDelimStrings[i].split(Delimiter)
        i = i+1

#=========================================================================

def groupByField (listOfDelimdStrings, fieldGroupedBy, sep='\t'): 
    """ WARNING!! listOfDelimdStrings will be destroyed!!
    Example:
    fieldGroupedBy = 0
    listOfDelimdStrings = ['a\t1','a\t2','b\t1']
    result-> [[['a','1'],['a','2']],[['b','1']] ]
    """

    # Convert listOfDelimdStrings to a list of lists with 'lowest' list
    # representing a list of original tab seperated fields
    listOfListsByTab = []
    while listOfDelimdStrings:
        #  Remove newLine and explode string on \t then append to listOfListsByTab
        listOfListsByTab.append(listOfDelimdStrings.pop(0).rstrip('\n').split(sep))

    exonList = []

    listOfExonLists = []

    while listOfListsByTab != 'end':

        # If there is a fresh and clean exonList:
        #     add the first coding region of the first/next gene to exonList
        if exonList == []:
            exonList.append(listOfListsByTab.pop(0))

            #  If that was the last entry, add it to listOfExonLists and set while loop up to end
            if listOfListsByTab == []:
                listOfExonLists.append(exonList)
                #print len(listOfExonLists), '\n'
                print exonList[0][fieldGroupedBy]
                exonList = []
                listOfListsByTab = 'end'


        # If the next BioMart record list matches the one(s) in exonList:
        #     add it to exonList
        elif listOfListsByTab[0][fieldGroupedBy] == exonList[0][fieldGroupedBy]:
            exonList.append(listOfListsByTab.pop(0))

            # Check to see if you just popped the last record:
            #    - export last exonList
            #    - cull exonList 
            #    - set listOfListsByTab to 'end' to stop the loop
            if listOfListsByTab == []:
                listOfExonLists.append(exonList)
                #print len(listOfExonLists), '\n'
                print exonList[0][fieldGroupedBy]
                exonList = []
                listOfListsByTab = 'end'

        # Otherwise append whole exonList to listOfExonLists and clean exonList for next record group
        else:
            listOfExonLists.append(exonList)
            #print len(listOfExonLists), '\n'
            print exonList[0][fieldGroupedBy]
            exonList = []

    print 'The groupByField function produced ',len(listOfExonLists),' groups.\n\n'

    return listOfExonLists


def groupByField_silent(listOfDelimdStrings, fieldGroupedBy, sep='\t'): 
    """ WARNING!! listOfDelimdStrings will be destroyed!!
    Example:
    fieldGroupedBy = 0
    listOfDelimdStrings = ['a\t1','a\t2','b\t1']
    result-> [[['a','1'],['a','2']],[['b','1']] ]
    """

    # Convert listOfDelimdStrings to a list of lists with 'lowest' list
    # representing a list of original tab seperated fields
    listOfListsByTab = []

    while listOfDelimdStrings != []:
        #  Remove newLine and explode sting on \t then append to listOfListsByTab
        listOfListsByTab.append(listOfDelimdStrings.pop(0).rstrip('\n').split(sep))
        lenListOfListsByTab = len(listOfListsByTab)
    exonList = []

    listOfExonLists = []
    lenListOfExonLists = len(listOfExonLists)
    while listOfListsByTab != 'end':

        # If there is a fresh and clean exonList:
        #     add the first coding region of the first/next gene to exonList
        if exonList == []:
            exonList.append(listOfListsByTab.pop(0))

            #  If that was the last entry, add it to listOfExonLists and set while loop up to end
            if listOfListsByTab == []:
                listOfExonLists.append(exonList)
                #print len(listOfExonLists), '\n'
                #print exonList[0][fieldGroupedBy]
                exonList = []
                listOfListsByTab = 'end'


        # If the next BioMart record list matches the one(s) in exonList:
        #     add it to exonList
        elif listOfListsByTab[0][fieldGroupedBy] == exonList[0][fieldGroupedBy]:
            exonList.append(listOfListsByTab.pop(0))

            # Check to see if you just popped the last record:
            #    - export last exonList
            #    - cull exonList 
            #    - set listOfListsByTab to 'end' to stop the loop
            if listOfListsByTab == []:
                listOfExonLists.append(exonList)
                #print len(listOfExonLists), '\n'
                print exonList[0][fieldGroupedBy]
                exonList = []
                listOfListsByTab = 'end'

        # Otherwise append whole exonList to listOfExonLists and clean exonList for next record group
        else:
            listOfExonLists.append(exonList)
            #print len(listOfExonLists), '\n'
            print exonList[0][fieldGroupedBy]
            exonList = []

    print 'The groupByField function produced ',len(listOfExonLists),' groups.\n\n'

    return listOfExonLists
