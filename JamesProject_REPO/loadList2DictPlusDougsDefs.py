#!/usr/bin/python
import math


#========================= User Defined Variables =========================
geneVectors     = 'geneVectors.txt'

listOfGeneNames = 'testCluster.txt'

#==========================================================================


#geneVectors = map(lambda line : line.strip().split('\t'), open(geneVectors, 'rU').readlines())

#cluster = map(lambda line : line.strip(), open(listOfGeneNames, 'rU').readlines())

## convert geneVectors to a_dict
#a_dict = {}

#for each in geneVectors:

	#listOfExpressionValues = each[1:]
	#for i in range(len(listOfExpressionValues)):
		#listOfExpressionValues[i] = float(listOfExpressionValues[i])
	
	#a_dict[each[0]] = listOfExpressionValues
	

a_dict = {'agap001': [1, 10, 4, 7], 'agap002': [2, 3, 4, 8], 'agap003': [1, 2, 5, 8], 
'agap004': [0, 2, 5, 6]}

cluster = ['agap001', 'agap003', 'agap004']	
	

cluster_vector_1 = []
for k in a_dict.keys():
	for i in cluster:
		if k == i:
			cluster_vector_1.append(a_dict[k])
cluster_vector_2 = []
for k in a_dict.keys():
	for i in cluster:
		if k == i:
			cluster_vector_2.append(a_dict[k])

def distance(list_1, list_2):
	sum_sq = 0 # sum of squares of the difference
	if len(list_1) != len(list_2):
		raise ValueError("Points have different dimensions")
	for i in range(len(list_1)):
		x, y = list_1[i], list_2[i] # extract the two ordinates
		sum_sq = sum_sq + (y-x)**2 # and add the square of the difference
		sqrt_sum_sq = math.sqrt(sum_sq)
	return sqrt_sum_sq

def PCC(vector_1, vector_2):
	sum_1 = 0
	aver_1 = 0
	for i in vector_1:
		sum_1 = sum_1 + i
		aver_1 = float(sum_1)/float(len(vector_1))
	sum_2 = 0
	aver_2 = 0
	for i in vector_2:
		sum_2 = sum_2 + i
		aver_2 = float(sum_2)/float(len(vector_2))
	diff_v1 = []
	for i in vector_1:
		diff = i - aver_1
		diff_v1.append(diff)
	diff_v2 = []
	for i in vector_2:
		diff = i - aver_2
		diff_v2.append(diff)
	multiplied = []
	for i in range(len(diff_v1)):
		x, y = diff_v1[i], diff_v2[i] # extract the two ordinates
		multiply = x * y
		multiplied.append(multiply)
	top_r = 0
	for i in multiplied:
		top_r = top_r + i
	sum_sqrt = []
	for i in diff_v1:
		squared = (i**2)
		sum_sqrt.append(squared)
	sum_sqrt_2 = []
	for i in diff_v2:
		squared = i**2
		sum_sqrt_2.append(squared)
	bottom_v1 = 0
	for i in sum_sqrt:
		bottom_v1 = bottom_v1 + i
		bottom_v2 = 0
		for j in sum_sqrt_2:
			bottom_v2 = bottom_v2 + j

	multi_bottom = (math.sqrt(bottom_v1*bottom_v2))

	r = top_r/multi_bottom

	return r

a = cluster_vector_1
b = cluster_vector_2

Euc_dists = []
for i in a:
	for j in b:
		temp = distance(i, j)
		Euc_dists.append(temp)

Pearsons = []
for i in a:
	for j in b:
		temp = PCC(i, j)
		Pearsons.append(temp)

print a, "\n", b, "\n", Euc_dists, Pearsons
