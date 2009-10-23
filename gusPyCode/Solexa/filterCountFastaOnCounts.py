#!/usr/bin/env python

import sys

#extract the reads with count equal or more than 500 from the file
#type: s_7_sequence_seqCount.fas

args = sys.argv[1:]

usage = "USAGE: python filterCountFastaOnCounts.py <inputFile> <outFile> <searchCount>"

assert len(args) == 3, usage

inFile = open(arg[0]), 'rU'

inFile  = open(args[0],'rU')
outFile = open(args[1], 'w')


while 1:
    line = inFile.readline()
    if not line: break
    # If line starts with '>' extract count and test it against user suppiled value (greater or equal)
    if line.startswith('>'):
        splitLine = line.strip().split('_')    # split ">0000000XX_N\n" on '_' after removing the \n
        # If count passes test add header
        # line AND seq line to outFile.
        if int(splitLine[1]) >= int(args[2]):  # must convert args[2] from text to int
            outFile.write(line)                # write 'line' to outFile
            outFile.write(inFile.readline())   # get next line (seqLine) and write it to outFile
            
inFile.close()
outFile.flush()
outFile.close()

