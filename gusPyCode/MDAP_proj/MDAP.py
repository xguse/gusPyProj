"""Usage: ???

MDAP_BASE is the base module for the Motif Discovery and Assesment Pipeline.
"""
#-------------------------
#        IMPORTS
#-------------------------
import sys
from optparse import OptionParser
from time import time


class MDAP():
    def __init__(self,argv):
	self.optParser = None	      # an OptionParser
	self.options   = None	      # the parsed options
	self.args      = None	      # the parsed positional args
	self.stdout    = sys.stdout   # the output file desc
	self.stderr    = sys.stderr   # the log file desc
	
	# set up the option parser's rules
	self.initOptParser()
	
	# send user's options through the grinder
	(self.options,self.args) = self.optParser.parse_args(argv)
	
	# validate user args supplied and further process any options that need it.
	self.validate()
	
	
	
	pass
    
    def initOptParser(self):
	usage = 'usage: %prog [options] fastaOfAllSeqs listOfLinkedSeqs nameForResults'
	self.optParser = OptionParser(usage=usage)
	# define options
	self.optParser.add_option('-k', '--kSize', dest='kmerSize', metavar='LEN', type='int', default=6, \
				  help='Sets motif length to use. [default=%default].')
	# NEW!!!
	self.optParser.add_option('-r', '--kRange', dest='kmerRange', metavar='LOW,HIGH', type='str', default=False, \
				  help='Sets motif range to use. [default=%default].')
	
	# --bg may be useless.  I can just use whether there is a 'memeBackground' provided or not.  AlignAce GC content is
	# fast enough to just ALWAYS calculate.
	self.optParser.add_option('--bg', dest='background', metavar='0/1', type='int', default=1, \
				  help="Compute each program's background requirement from data. [default=%default]") 
	
	self.optParser.add_option('--mBg', dest='memeBackground', metavar='FILE', type='string', default=None, \
				  help="File with an n-order Markov background model for Meme. [default=%default]")
	
	self.optParser.add_option('--aiter', dest='ace_iter', metavar='NUM', type='int', default=3, \
				  help='Iterations of AlignAce to run. [default=%default]')
	
	self.optParser.add_option('-a', '--ace', dest='alignAce', action='store_true', default=False, \
				  help='Include AlignAce results. [default=%default]')
	
	self.optParser.add_option('-m', '--meme', dest='meme', action='store_true', default=False, \
				  help='Include MEME results. [default=%default]')
	
	self.optParser.add_option('-s', '--mdscan', dest='mdScan', action='store_true', default=False, \
				  help='Include MDScan results. [default=%default]')
	
	self.optParser.add_option('-t', '--tamo', dest='TAMO', action='store_true', default=False, \
				  help='Include TAMO results. [default=%default]')
	
	# to add: markovBgFile, 
	
	
	
    
    def validate(self):
	"""Validate user options and returns helpful error messages."""
	# NOTES: 
	# ___ Convert options['kmerRange'] from str<LOW,HIGH> to list of ints.
	# ___ If True in list<appSwitches>(must create): only use those switches with True, else: use all apps
	pass
    
    def go(self):
	"""Main execution function.  Initializes data sets and calls the separate discovery programs.
	Uses parsed options to determine what to do with the output."""
	pass


if __name__ == "__main__":
	MDAP(sys.argv).go()
