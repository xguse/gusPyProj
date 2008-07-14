#!/usr/bin/python

import sys


list_1 = open('cluster_1.txt', 'r').readlines()
list_2 = open('cluster_2.txt', 'r').readlines()




i = 0
while i < len(list_1):
    list_1[i] = list_1[i].rstrip('\n')
    i = i + 1

j = 0
while j < len(list_2):
    list_2[j] = list_2[j].rstrip('\n')
    j = j + 1

#print set_1


list1_Only = []
list2_Only = []
union = []


for list1_Probe in list_1:
    if list1_Probe not in list_2:
        list1_Only.append(list1_Probe)
    else:
        union.append(list1_Probe)

for list2_Probe in list_2:
    if list2_Probe not in list_1:
        list2_Only.append(list2_Probe)
        
print 'yay!'
        
