#!/usr/bin/python

import sys
import string
motifs = open(sys.argv[1], 'r').readlines()

cleaned = map(string.strip, motifs)
#print cleaned

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

compl_iupacdict = {'A':'T',
	'C':'G',
	'G':'C',
	'T':'A',
	'M':'[TG]',
	'R':'[TC]',
	'W':'[TA]',
	'S':'[GC]',
	'Y':'[GA]',
	'K':'[CA]',
	'V':'[TGC]',
	'H':'[TGA]',
	'D':'[TCA]',
	'B':'[GCA]',
	'X':'[ACGT]',
	'N':'[ACGT]'}

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

def reverse(text):
	return ''.join([text[i] for i in range(len(text)-1,-1,-1)])


for i in cleaned:
	s_listings = iupac_transl(i, iupacdict)
	c_listings = iupac_transl(i, compl_iupacdict)
	co_listings = reverse(c_listings)
	print i, "\t", s_listings, "\t", c_listings, "\t", co_listings