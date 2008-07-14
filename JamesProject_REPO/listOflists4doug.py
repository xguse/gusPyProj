#========================= User Defined Variables =========================
pathToFile      = 'dougListTest.txt'  #'path/to/file' 

#==========================================================================


#  this gives you a list of the lines in the file with the \n's stripped
inFile = map(lambda line: line.strip(), open(pathToFile, 'rU').readlines())

#  now we need to split all of the list items which are strings at the moment into lists
for index in range(0,len(inFile)):
    #  this replaces the item that is a string with the list version for each list index
    inFile[index] = inFile[index].split('\t')
    
    #  now we float-ify each number in the list since they are strings right now
    for number in range(0,len(inFile[index])):
        #  this notation is how you access items that are in the 2dary lists:  list[i][j]
        inFile[index][number] = float(inFile[index][number])

#  now we do a REAL copy of the list in the form of  l2 = l1[:]
#  this slices the first list's indexs to the new list, but nothing is destroyed
#  so you end up with an index by index REAL copy of l1
inFileCopy = inFile[:]  



