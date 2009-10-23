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
	'N':'N'}

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



	
def iupac2regex(motif):

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

def groupByField (listOfTabbedStrings, fieldGroupedBy): 
	""" WARNING!! listOfTabbedStrings will be destroyed!!
    Example:
    fieldGroupedBy = 0
    listOfTabbedStrings = ['a\t1','a\t2','b\t1']
    result-> [[['a','1'],['a','2']],[['b','1']] ]
    """

	# Convert listOfTabbedStrings to a list of lists with 'lowest' list
	# representing a list of original tab seperated fields
	listOfListsByTab = []
	while listOfTabbedStrings != []:
		#  Remove newLine and explode sting on \t then append to listOfListsByTab
		listOfListsByTab.append(listOfTabbedStrings.pop(0).rstrip('\n').split('\t'))

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

    
def groupByField_silent(listOfTabbedStrings, fieldGroupedBy): 
	""" WARNING!! listOfTabbedStrings will be destroyed!!
    Example:
    fieldGroupedBy = 0
    listOfTabbedStrings = ['a\t1','a\t2','b\t1']
    result-> [[['a','1'],['a','2']],[['b','1']] ]
    """

	# Convert listOfTabbedStrings to a list of lists with 'lowest' list
	# representing a list of original tab seperated fields
	listOfListsByTab = []
	
	while listOfTabbedStrings != []:
		#  Remove newLine and explode sting on \t then append to listOfListsByTab
		listOfListsByTab.append(listOfTabbedStrings.pop(0).rstrip('\n').split('\t'))
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
			#print exonList[0][fieldGroupedBy]
			exonList = []

	print 'The groupByField function produced ',len(listOfExonLists),' groups.\n\n'

	return listOfExonLists
