#! /usr/bin/python

class Fasta:
	def __init__(self, name, sequence):
		self.name = name  # this will store the sequence name
		self.sequence = sequence # this will store the sequence


#this function will receive the list with the file
#contents, create instances of the Fasta class as
#it scans the list, putting the sequence name on the
#first attribute and the sequence itself on the second
#attribute

def read_fasta(file):
	items = []
	index = 0
	for line in file:
		if line.startswith(">"): # we check to see if the line starts with a > sign
			#if so and our counter is large than 1 we add the created class instance
			#to our list a counter larger than 1 means we are reading from sequences
			#2 and above
			if index >= 1:
				items.append(aninstance)
			index+=1
			#we add the line contents to a string
			name = line[:-1]
			#and initialize the string to store the sequence
			seq = ''
			#this creates a class instance and we add the attributes
			#which are the strings name and seq
			aninstance = Fasta(name, seq)
		else:
			#the line does not start with > so it has to be
			#a sequence line, so we increment the string and 
			#add it to the created instance
			seq += line[:-1]
			aninstance = Fasta(name, seq)

	#the loop before reads everything but the penultimate
	#sequence is added at the end, so we need to add it
	#after the loop ends
	items.append(aninstance)
	#a list with all read sequences is returned
	return items

#fastafile = open("test.fas", "r").readlines()
#mysequences = read_fasta(fastafile)

#print mysequences

#for i in mysequences:
	#print i.name
#    print i.sequence


def format_output(sequence, length):
	temp = []
	for j in range(0,len(sequence),length):
		temp.append(sequence[j:j+length])
	return '\n'.join(temp)