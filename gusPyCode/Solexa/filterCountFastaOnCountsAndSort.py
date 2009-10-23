#!/usr/bin/env python

import sys

args = sys.argv[1:]

usage = """\n\nUSAGE: python sortCountFasta.py inFile outFile "sortDefinition"

\tsortDefinition examples: (NOTE: you MUST use 'x' and bound the definition in "")
\t5<x<=300
\tx>=2000"""

assert len(args) == 3, usage
print 'Sort Instructions = %s' % (args[-1])

inPath  = args[0]
outPath = args[1]
inFile  = open(inPath,'rU')
outFile = open(outPath,'w')

data = []

print 'Looping over file...'
while 1:
    line = inFile.readline()
    if not line: break
    # If line starts with '>' take it AND the next line and create a single list representing all the parts
    #     example: [id, count, seq]
    # Store that list in the list named 'data'
    tempList = [] # Create temperary list to store fasta data
    if line.startswith('>'):
        x = int(line.strip('\n').split('_')[-1])                        # get count for seq
        #print x
        if eval(args[2]):                                               # use sortDefinition to decide if the fasta should be kept
            #print 'kept '+str(x)
            tempList.extend(line.strip('\n').lstrip('>').split('_'))    # split ">0000000XX_N\n" on '_' after removing the \n and the > then add to tempList (lstrip() strips on the left of teh str)
            tempList.append(inFile.readline().strip())                  # append seq to tempList
            data.append(tempList)                                       # append a copy of tempList to data

# we're done with the infile for now so close it 
inFile.close()
        
# Sort data using data.sort() but with a fancy definition that tells
# data to sort based on the second index of each list in data (the count item)
data.sort(key=lambda x: int(x[1]))  # lamda is a special tool to define a simple subroutine on the fly DONT WORRY ABOUT UNDERSTANDING IT. 
data.reverse()                      # sort() sorts ascending by default.  So We reverse it to get the biggest atthe top. 
            


for line in data:
    # formats the data back into the same fasta form as
    # we found it in and writes it to the file
    outFile.write('>%s_%s\n%s\n' % (line[0],line[1],line[2]))  
    
outFile.flush()
outFile.close()

print 'I have finished.'
