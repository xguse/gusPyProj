"""This is a short script that runs MDAP's AceWrap wapper.  Can be used to show that the
system is correctly set up to run AlignAce."""

from MD_wrappers import AceWrap
import sys

# AceWrap -> def __init__(self, optionsObj, posArgs):

optionsObj = {'kmerSize':8,
              'kmerRange':'6,9',
              'background':1,
              'ace_iter':10,} 

posArgs    = ['/Users/biggus/Documents/James/Data/2KB/2kb_Sequence/2kb_Anopheles/2KBupTSS_goodAffyAGAPsFastasOUT.masked.nr.fas',
              '/Users/biggus/Documents/James/Data/DoubleBloodMeal/DBMclusterGenes.txt']

outFile    = '/Users/biggus/Documents/James/Data/DoubleBloodMeal/DBMgenesAAceAnalysis.8mers.txt'

aWrap = AceWrap(optionsObj,posArgs)
aWrap.go()

#for l in tWrap.toFile:
    #print l.strip()

print 'num of motifs = %s' % (len(aWrap.output.results))


# trickiness to capture print to variable
sys.stdout = open(outFile, 'w')
count = 1
for m in aWrap.output.results:
    print 'AlignAce Results for motif %s' % (count)
    count+=1
    print m
    print'\n>%s\n%s\n' % (m.oneletter,m.oneletter)
    print'Probability Matrix:'
    m._print_p()
    print '\nText Logo:'
    m.printlogo()
    print '\n\n***************************************************\n\n'
    
print 'Raw AlignAce Results from each iteration:\n%s' % (''.join(aWrap.output.lines))

# reset stdout
sys.stdout = sys.__stdout__

print 'All done.'

None

