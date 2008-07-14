from sets import Set
from time import time

#========================= User Defined Variables =========================
pValsFile      = '/Users/biggus/Documents/MBGB/Rotations/James/Data/2KB/2kb_ClusterHyperGeo/comboHyperGeo/timeCourseMotifCombos.pvals.txt'

comboName      = '12_time-course*3h_UP'

outFile        = '/Users/biggus/Documents/MBGB/Rotations/James/Data/2KB/2kb_ClusterHyperGeo/comboHyperGeo/%s.pvals.txt' % (comboName)

#===============================================

pValsFile = map(lambda line: line.strip(), open(pValsFile, 'rU').readlines())

listOfPvalsForCluster = []

for each in pValsFile:
    firstField = each.split('\t')[0].split(';')
    if firstField[1] == comboName:
        listOfPvalsForCluster.append(each+'\n')
    
outFile = open(outFile, 'w').writelines(listOfPvalsForCluster)

print 'cluster/combos = %s lines.' % (str(len(listOfPvalsForCluster)))
