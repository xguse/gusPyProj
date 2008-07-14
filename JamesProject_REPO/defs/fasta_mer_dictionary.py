#!/usr/bin/python/

import fasta
import sys

fastafile = open("test.fas", "r").readlines()
my_sequences = fasta.read_fasta(fastafile)

#test_list = ['ATATAG', 'TATA', 'GGGTGA']

#to make a dictionary of all possible hexamers
all_6mer = {}
base = ['A','C','G','T']
for C1 in base:
	for C2 in base:
		for C3 in base:
			for C4 in base:
				for C5 in base:
					for C6 in base:
						all_6mer[''.join([C1,C2,C3,C4,C5,C6])] = 0

length = len(all_6mer)

keys = all_6mer.keys()

for i in my_sequences:
	ind_sequence = i.sequence
	print i.name
	for keys in all_6mer:
		match = {}
		if keys in ind_sequence:
			print keys
		#else:
				#print keys,"\t", 0


#for i in my_sequences:
#    ind_sequence = i.sequence
#    print i.name
#    for i in test_list:
#        if i in ind_sequence:
#           print i, 1
#        else:
#             print i, 0






#for keys in all_6mer:
#    if keys in seq_test:
#       print keys
