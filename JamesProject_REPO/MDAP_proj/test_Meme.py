"""This is a short script that runs MDAP's AceWrap wapper.  Can be used to show that the
system is correctly set up to run AlignAce."""

from MD_wrappers import MemeWrap
import sys

# AceWrap -> def __init__(self, optionsObj, posArgs):

optionsObj = {'kmerSize':8,
              'kmerRange':'6,8',
              'background':1,
              'memeBackground':None} 

posArgs    = ['/Users/biggus/Documents/James/Data/2KB/2kb_Sequence/2kb_Anopheles/2KBupTSS_goodAffyAGAPsFastasOUT.masked.nr.fas',
              '/Users/biggus/Documents/James/Data/DoubleBloodMeal/DBMclusterGenes.testMeme.txt']

outFile    = '/Users/biggus/Documents/James/Data/DoubleBloodMeal/DBMgenesMemeAnalysis.6-7mers.testMeme.txt'

mWrap = MemeWrap(optionsObj,posArgs)
mWrap.go()

#for l in tWrap.toFile:
    #print l.strip()

print 'num of motifs = %s' % (len(mWrap.output.motifs))


## trickiness to capture print to variable
#sys.stdout = open(outFile, 'w')
#count = 1
#for m in mWrap.output.results:
    #print 'Meme Results for motif %s' % (count)
    #count+=1
    #print m
    #print'\n>%s\n%s\n' % (m.oneletter,m.oneletter)
    #print'Probability Matrix:'
    #m._print_p()
    #print '\nText Logo:'
    #m.printlogo()
    #print '\n\n***************************************************\n\n'
    
#print 'Raw AlignAce Results from each iteration:\n%s' % (''.join(aWrap.output.lines))

## reset stdout
#sys.stdout = sys.__stdout__

print 'All done.'

None

